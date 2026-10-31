"""PlateSim, plate events: splitting, merging, detached pieces and welds (B-b), rifting of plates and subduction
initiation."""
import numpy as np

from planet_grid import fields

from planet_tectonics.graphops import components, nearest_label
from planet_tectonics.params import COLLISION


class EventsMixin:
    """PlateSim, plate events: splitting, merging, detached pieces and welds (B-b), rifting of plates and subduction
    initiation."""

    def split(self, pl, part_world, rel_speed_cm_yr, kind):
        """Split plate `pl`: the world nodes in part_world (owned by pl) become a new plate in the same frame (no
        resampling). The two parts get Euler vectors differing by rel_speed (positive: diverge, negative: converge),
        about the pole 90 degrees from both parts' centroids."""
        g = self.g
        own = self.owner == pl.id
        part = part_world & own
        rest = own & ~part
        if part.sum() < 4 or rest.sum() < 4:
            return None
        wn = fields.nearest_nodes(g, g.xyz@pl.Q.T)
        lab = part[wn] & pl.mask
        new = self._new_plate(pl.axis, pl.speed)
        new.Q = pl.Q.copy()
        for name in ('mask', 'cont', 'thick', 'born', 'orog', 'under', 'sed', 'slab'):
            arr = getattr(pl, name)
            na = getattr(new, name)
            na[lab] = arr[lab]
        pl.mask[lab] = False; pl.cont[lab] = False
        c1 = (g.area[rest, None]*g.xyz[rest]).sum(0); c1 /= np.linalg.norm(c1)
        c2 = (g.area[part, None]*g.xyz[part]).sum(0); c2 /= np.linalg.norm(c2)
        d = np.cross(c1, c2)
        if np.linalg.norm(d) < 1e-9:
            d = self.rng.normal(size=3)
        d = d/np.linalg.norm(d)*rel_speed_cm_yr
        base = self._omega_vec(pl)
        self._set_omega(pl, base-d/2)
        self._set_omega(new, base+d/2)
        self.owner = np.where(part, new.id, self.owner)
        self.events.append({'t': self.t, 'kind': kind, 'plate': int(pl.id), 'new_plate': int(new.id),
                            'area_mkm2': float(g.area[part].sum()/1e12)})
        return new


    def _merge_events(self):
        g, p = self.g, self.p
        coll = np.flatnonzero(self.btype == COLLISION)
        contact = {}
        if coll.size:
            nb = g.nbr[coll]
            ok = nb >= 0
            a = np.repeat(self.owner[coll][:, None], nb.shape[1], 1)
            b = self.owner[np.where(ok, nb, 0)]
            m = ok & (a != b) & (self.btype[np.where(ok, nb, 0)] == COLLISION)
            pairs, counts = np.unique(np.stack((np.minimum(a[m], b[m]), np.maximum(a[m], b[m])), 1), axis=0, return_counts=True)
            h = np.sqrt(g.area.mean())/1e3
            for (x, y), c in zip(pairs, counts):
                contact[(int(x), int(y))] = c*h/2          # each contact appears from both sides
        prev = getattr(self, 'contact_myr', {})
        self.contact_myr = {k: prev.get(k, 0.)+p['dt_myr'] for k, L in contact.items() if L >= p['merge_contact_km']}
        for (x, y), myr in list(self.contact_myr.items()):
            if myr < p['merge_after_myr']:
                continue
            pa = [q for q in self.plates if q.id == x]; pb = [q for q in self.plates if q.id == y]
            if not pa or not pb:
                continue
            ax, ay = g.area[self.owner == x].sum(), g.area[self.owner == y].sum()
            big, small = (pa[0], pb[0]) if ax >= ay else (pb[0], pa[0])
            if min(ax, ay) > p['merge_max_ratio']*max(ax, ay) or (ax+ay)/1e12 > p['merge_max_area_mkm2']:
                continue
            self.merge(big, small)
            self.contact_myr = {k: v for k, v in self.contact_myr.items() if small.id not in k}


    def merge(self, big, small):
        """Weld `small` onto `big`. Every node of big's raster whose world position lies in small's territory reads
        small's crust there (backward lookup, so the region is filled without holes); where big already holds crust,
        small's continental volume is added to big's thickness. Euler vectors are averaged by area."""
        g = self.g
        v0 = self.continental_volume()
        vs0 = self.continental_sediment_volume()
        wk = fields.nearest_nodes(g, g.xyz@big.Q.T)                  # big node -> world node
        take = self.owner[wk] == small.id
        k = np.flatnonzero(take)
        j = fields.nearest_nodes(g, (g.xyz[k]@big.Q.T)@small.Q)      # -> small's raster
        ok = small.mask[j]
        k, j = k[ok], j[ok]
        free = ~big.mask[k]
        for name in ('cont', 'thick', 'born', 'orog', 'sed', 'slab'):
            getattr(big, name)[k[free]] = getattr(small, name)[j[free]]
        big.under[k[free]] = 0
        big.mask[k[free]] = True
        # the old boundary is now inside one plate: every big node whose world position the welded plate owns must
        # hold crust, else later sampling finds cracks along the seam; fill them from the nearest set node
        need = ((self.owner[wk] == big.id) | take) & ~big.mask
        if need.any():
            src_lab = np.where(big.mask, np.arange(g.size), -1).astype(np.int64)
            near, _ = nearest_label(src_lab, g.nbr, g.dist)
            src = near[need]
            big.mask[need] = True
            for name in ('cont', 'thick', 'born', 'orog', 'sed', 'slab'):
                arr = getattr(big, name); arr[need] = arr[src]
        both = ~free & small.cont[j] & big.cont[k]
        self.area_ledger['merge'] = self.area_ledger.get('merge', 0.)-float(g.area[k[both]].sum())
        np.add.at(big.thick, k[both], small.thick[j[both]])
        if self.p['continental_sediment']:
            np.add.at(big.sed, k[both], small.sed[j[both]])
        a_big, a_small = g.area[self.owner == big.id].sum(), g.area[self.owner == small.id].sum()
        self._set_omega(big, (self._omega_vec(big)*a_big+self._omega_vec(small)*a_small)/(a_big+a_small))
        self.owner = np.where(self.owner == small.id, big.id, self.owner)
        self.plates = [q for q in self.plates if q is not small]
        self.ledger['merge'] = self.ledger.get('merge', 0.)+self.continental_volume()-v0
        self.sed_ledger['merge'] += self.continental_sediment_volume()-vs0
        self.events.append({'t': self.t, 'kind': 'merge', 'plate': int(big.id), 'merged': int(small.id),
                            'area_mkm2': float(a_small/1e12)})


    # ------------------------------------------------------------------ plate pieces (B-b)
    def plate_pieces(self):
        """Every plate's world territory as connected pieces: (plate, component labels, piece areas m^2, index of
        the largest piece) for each plate that owns nodes."""
        g = self.g
        out = []
        for pl in self.plates:
            own = self.owner == pl.id
            if not own.any():
                continue
            lab = components(own, g.nbr)
            sz = np.bincount(lab[lab >= 0], weights=g.area[lab >= 0])
            out.append((pl, lab, sz, int(np.argmax(sz))))
        return out


    def detached_area(self):
        """World area (m^2) of the plate pieces other than each plate's largest."""
        return float(sum(sz.sum()-sz[m] for _, _, sz, m in self.plate_pieces()))


    def _separate_pieces(self):
        """piece_rule 'plates' (B-b): each detached piece with an interior node (all its neighbours on the same
        plate) becomes its own plate if it covers at least piece_min_area_mkm2, else joins the plate around most of
        its edge (_weld)."""
        g, p = self.g, self.p
        floor = p['piece_min_area_mkm2']*1e12
        nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
        st = self.piece_stats = getattr(self, 'piece_stats', {'plates_made': 0, 'welds': 0, 'welded_m2': 0.,
                                                             'crumbs_m2_steps': 0.})
        small = []
        for pl, lab, sz, main in self.plate_pieces():
            if sz.size < 2:
                continue
            ids_all = np.flatnonzero(lab >= 0)
            o = ids_all[np.argsort(lab[ids_all], kind='stable')]
            cuts = np.searchsorted(lab[o], np.arange(sz.size+1))
            for c in range(sz.size):
                if c == main:
                    continue
                ids = o[cuts[c]:cuts[c+1]]
                if not (self.owner[nb[ids]] == pl.id).all(1).any():
                    st['crumbs_m2_steps'] += float(g.area[ids].sum())
                    continue
                if sz[c] >= floor:
                    part = np.zeros(g.size, bool); part[ids] = True
                    new = self.split(pl, part, 0., 'detach')
                    if new is not None:
                        self._wn[new.id] = self._wn[pl.id]          # same frame: same raster -> world map
                        st['plates_made'] += 1
                        continue
                small.append((pl.id, ids))
        groups = {}
        for sid, ids in small:                       # after the splits: the surrounding plate may be a new one
            nn = g.nbr[ids]
            o = self.owner[np.where(nn >= 0, nn, 0)]
            ok = (nn >= 0) & (o != sid)
            if not ok.any():
                continue
            u, cnt = np.unique(o[ok], return_counts=True)
            groups.setdefault((sid, int(u[np.argmax(cnt)])), []).append(ids)
        byid = {q.id: q for q in self.plates}
        for (sid, did), parts in groups.items():
            if sid in byid and did in byid:
                ids = np.concatenate(parts)
                self._weld(byid[sid], byid[did], ids)
                st['welds'] += 1
                st['welded_m2'] += float(g.area[ids].sum())


    def _weld(self, src, dst, ids):
        """Move the world nodes ids (a piece of src) into dst's raster, as merge does for a whole plate: dst's raster
        nodes whose world position is on the piece read src's crust there (backward lookup); where dst already holds
        crust under it, src's crust replaces dst's oceanic crust (subducted) and stacks on its continental crust
        (collision, as merge); dst keeps its motion. Volume changes go to the 'pieces' ledgers."""
        g = self.g
        a = g.area
        part = np.zeros(g.size, bool); part[ids] = True
        wd, ws = self._wn[dst.id], self._wn[src.id]
        k = np.flatnonzero(part[wd])                     # dst raster nodes on the piece
        ro = np.flatnonzero(src.mask & part[ws])         # src raster nodes on the piece

        def vol(pl, r, f):
            return float((f[r].astype(np.float64)*a[r]*(pl.mask[r] & pl.cont[r])).sum())
        names = ('cont', 'thick', 'born', 'orog', 'sed', 'slab')
        v0 = vol(src, ro, src.thick)+vol(dst, k, dst.thick)
        s0 = vol(src, ro, src.sed)+vol(dst, k, dst.sed)
        if k.size:
            j = fields.nearest_nodes(g, (g.xyz[k]@dst.Q.T)@src.Q)
            ok = src.mask[j]
            k1, j1 = k[ok], j[ok]
            free = ~dst.mask[k1]
            take = free | ~dst.cont[k1]
            stack = ~take & src.cont[j1]
            kt, jt = k1[take], j1[take]
            for name in names:
                getattr(dst, name)[kt] = getattr(src, name)[jt]
            dst.mask[kt] = True
            dst.under[kt] = 0
            ks, js = k1[stack], j1[stack]
            np.add.at(dst.thick, ks, src.thick[js])
            if self.p['continental_sediment']:
                np.add.at(dst.sed, ks, src.sed[js])
            for _ in range(3):                           # cracks: copy from a set neighbour (as merge's seam fill)
                need = k[~dst.mask[k]]
                if not need.size:
                    break
                nn = g.nbr[need]
                m = (nn >= 0) & dst.mask[np.where(nn >= 0, nn, 0)]
                has = m.any(1)
                if not has.any():
                    break
                sn = np.take_along_axis(nn, np.argmax(m, 1)[:, None], 1)[:, 0]
                nd, sn = need[has], sn[has]
                for name in names:
                    arr = getattr(dst, name); arr[nd] = arr[sn]
                dst.mask[nd] = True
                dst.under[nd] = 0
        src.mask[ro] = False
        src.cont[ro] = False
        self.ledger['pieces'] = self.ledger.get('pieces', 0.)+vol(src, ro, src.thick)+vol(dst, k, dst.thick)-v0
        self.sed_ledger['pieces'] = self.sed_ledger.get('pieces', 0.)+vol(src, ro, src.sed)+vol(dst, k, dst.sed)-s0
        self.owner = np.where(part, dst.id, self.owner)


    def _split_events(self):
        g, p, rng, w = self.g, self.p, self.rng, self.world
        dt = p['dt_myr']
        A = g.area
        for pl in list(self.plates):
            own = self.owner == pl.id
            area = A[own].sum()/1e12
            if area < p['split_min_area_mkm2']:
                continue
            cf = A[own & w['cont']].sum()/A[own].sum()
            rate = p['split_rate_per_myr']*area/100.*(2. if cf > .5 else 1.)
            if rng.uniform() < 1-np.exp(-rate*dt):
                self._rift_plate(pl, own)
                continue
            old = own & ~w['cont'] & (self.t-w['born'] > p['init_age_myr'])
            old_area = A[old].sum()/1e12
            irate = p['init_rate_per_myr']*min(5., old_area/p['init_area_ref_mkm2'])
            if old.sum() >= 4 and rng.uniform() < 1-np.exp(-irate*dt):
                self._initiate_subduction(pl, own, old)


    def _rift_plate(self, pl, own):
        """Two seeds far apart inside the plate; the split follows a wiggled Voronoi line that prefers old sutures."""
        g, p, rng, w = self.g, self.p, self.rng, self.world
        ids = np.flatnonzero(own)
        a = ids[rng.integers(ids.size)]
        dist_a = np.arccos(np.clip(g.xyz[ids]@g.xyz[a], -1, 1))
        b = ids[rng.choice(np.argsort(dist_a)[-max(1, ids.size//10):])]
        weak = np.where(np.isnan(w['orog']), 0., p['suture_weakness'])
        wig = 1+.3*rng.standard_normal(g.size)
        cost = np.clip(.5*(wig[:, None]+wig[np.clip(g.nbr, 0, None)]), .3, None)*(1-.5*(weak[:, None]+weak[np.clip(g.nbr, 0, None)]))
        cost = np.where(own[:, None] & own[np.clip(g.nbr, 0, None)], g.dist*cost, 1e30)
        lab0 = np.full(g.size, -1, np.int64); lab0[a] = 0; lab0[b] = 1
        lab, _ = nearest_label(lab0, g.nbr, cost)
        self.split(pl, lab == 1, p['split_speed_cm_yr'], 'rift')


    def _initiate_subduction(self, pl, own, old):
        """Old ocean breaks away from the rest of its plate and converges on it."""
        g, p, w = self.g, self.p, self.world
        ocean = own & ~w['cont']
        lab = components(ocean, g.nbr)
        target = np.bincount(lab[old & (lab >= 0)]).argmax()
        part = lab == target
        if (own & ~part).sum() == 0:                     # all-ocean plate: split through the old region instead
            c = (g.area[old, None]*g.xyz[old]).sum(0); c /= np.linalg.norm(c)
            n = np.cross(c, self.rng.normal(size=3)); n /= np.linalg.norm(n)
            part = own & (g.xyz@n > 0)
        new = self.split(pl, part, -p['init_speed_cm_yr'], 'subduction initiation')
        if new is not None:                              # spontaneous initiation: the old slab sinks on its own (Stern 2004)
            self.slab_seed_pairs = getattr(self, 'slab_seed_pairs', set()) | {(min(pl.id, new.id), max(pl.id, new.id))}
