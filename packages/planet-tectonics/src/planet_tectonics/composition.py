"""PlateSim, composition: which plate owns each world node (claims, scores, ocean-ocean polarity), the world fields,
and the plate-raster update after a step (hidden and consumed crust, collision volume, increments, ledgers)."""
import numpy as np

from planet_grid import fields

from planet_tectonics.crust_kernels import collision_kernel
from planet_tectonics.graphops import nearest_label, plate_holes
from planet_tectonics.lithosphere import thickened_therm


class CompositionMixin:
    """PlateSim, composition: which plate owns each world node (claims, scores, ocean-ocean polarity), the world fields,
    and the plate-raster update after a step (hidden and consumed crust, collision volume, increments, ledgers)."""

    def _label_boundary(self, mask, lab):
        g = self.g
        nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
        return mask & ((mask[nb] & (lab[nb] != lab[:, None])).any(1))


    def _other_label(self, mask, lab):
        """For each node, the label of the nearest node of a different label within the mask."""
        g = self.g
        nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
        diff = mask[nb] & (lab[nb] != lab[:, None]) & mask[:, None]
        src = np.where(diff.any(1), np.take_along_axis(lab[nb], np.argmax(diff, 1)[:, None], 1)[:, 0], -1)
        out, _ = nearest_label(np.where(mask, src, -1), g.nbr, g.dist)
        return out


    # ------------------------------------------------------------------ composition
    def _claims(self):
        g = self.g
        K = len(self.plates)
        idx = np.empty((K, g.size), np.int64)
        claimed = np.empty((K, g.size), bool)
        for k, pl in enumerate(self.plates):
            idx[k] = fields.nearest_nodes(g, g.xyz@pl.Q)          # Q^T p, as row vectors
            claimed[k] = pl.mask[idx[k]]
        return idx, claimed


    def _score(self, cont, age, cont_share):
        """Overlap priority: continental crust (bigger continent first) beats ocean; younger ocean beats older."""
        return np.where(cont, 1e6+1e5*cont_share, -age)


    def _ocean_polarity(self, claimed, cont_k, age_k, nclaim, top):
        """Ocean-ocean overlaps: one subduction polarity per pair of plates, not per node. The plate whose crust along
        the contact is younger (lighter) overrides; the polarity is remembered while the pair stays in contact, as
        real subduction zones keep their polarity."""
        ocean_only = (nclaim >= 2) & ~(claimed & cont_k).any(0)
        if not ocean_only.any():
            self.polarity_now = {}
            return top
        pid = np.array([pl.id for pl in self.plates])
        pol = getattr(self, 'polarity', {})
        now = self.polarity_now = {}                     # this step's pairs (B-d2, trench_polarity 'pair'): (lo id, hi id) -> over
        cols = np.flatnonzero(ocean_only)
        cl = claimed[:, cols]
        K = len(self.plates)
        wins = np.zeros((K, cols.size), np.int32)
        for a in range(K):
            for b in range(a+1, K):
                both = cl[a] & cl[b]
                if not both.any():
                    continue
                key = (pid[a], pid[b])
                if key not in pol:
                    ya = age_k[a, cols[both]].mean(); yb = age_k[b, cols[both]].mean()
                    pol[key] = pid[a] if ya <= yb else pid[b]
                over = a if pol[key] == pid[a] else b
                now[(min(pid[a], pid[b]), max(pid[a], pid[b]))] = pol[key]
                wins[over] += both
        self.polarity = pol
        best = np.argmax(np.where(cl, wins, -1), 0)
        top = top.copy()
        top[cols] = best
        return top


    def compose_from_plates(self):
        """World fields from the plate rasters without moving anything (start, and after events)."""
        idx, claimed = self._claims()
        self._world_from(idx, claimed, record=False)


    def _world_from(self, idx, claimed, record=True):
        g, p = self.g, self.p
        K = len(self.plates)
        cols = np.arange(g.size)
        cont_k = np.stack([pl.cont[idx[k]] for k, pl in enumerate(self.plates)])
        age_k = np.stack([self.t-pl.born[idx[k]] for k, pl in enumerate(self.plates)])
        share = np.array([(g.area*(pl.mask & pl.cont)).sum() for pl in self.plates])/g.area.sum()
        sc = np.where(claimed, self._score(cont_k, age_k, share[:, None]), -np.inf)
        nclaim = claimed.sum(0)
        top = np.argmax(sc, 0)
        top = self._ocean_polarity(claimed, cont_k, age_k, nclaim, top)
        gap = nclaim == 0
        # gaps join the nearest plate (ridge); label by slot
        lab = np.where(gap, -1, top)
        lab, _ = nearest_label(lab.astype(np.int64), g.nbr, g.dist)
        slot = lab
        pid = np.array([pl.id for pl in self.plates])
        self.owner_slot = slot
        self.owner = pid[slot]
        if p['compose_map'] == 'matched':
            self._match_reads(idx, slot, gap)
            cont_k = np.stack([pl.cont[idx[k]] for k, pl in enumerate(self.plates)])
        elif p['compose_map'] != 'nearest':
            raise ValueError(f"compose_map {p['compose_map']!r}: 'nearest' or 'matched'")
        wcont = np.where(gap, False, cont_k[slot, cols])
        wthick = np.where(gap, p['ocean_thick_m'], np.stack([pl.thick[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        wborn = np.where(gap, self.t, np.stack([pl.born[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        worog = np.where(gap, np.nan, np.stack([pl.orog[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        wsed = np.where(gap, 0., np.stack([pl.sed[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        wtherm = np.where(gap, np.nan, np.stack([pl.therm[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        wrock = np.where(gap, 0., np.stack([pl.rock[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        if p['remap_rule'] in ('tent', 'overlap'):        # pass 108: interpolated / overlap read (see PARAMS)
            wthick, wsed, wrock, wborn, wtherm = self._tent_gather(idx, slot, gap, wthick, wsed, wrock, wborn, wtherm)
        elif p['remap_rule'] != 'nearest':
            raise ValueError(f"remap_rule {p['remap_rule']!r}: 'nearest', 'tent' or 'overlap'")
        wrock = np.clip(wrock, 0., np.maximum(wsed, 0.))         # sediment_rock: never more rock than sediment
        self.new_crust = gap
        self.world = {'cont': wcont, 'thick': wthick, 'born': wborn, 'orog': worog, 'sed': wsed, 'therm': wtherm, 'rock': wrock}
        self.therm_changed = np.zeros(g.size, bool)      # C: world nodes whose thermal age changed this step


    def _tent_weights(self, pl, ws, j0):
        """remap_rule 'tent' (pass 108, see PARAMS): for world nodes ws read from plate pl at nearest plate nodes j0, the
        candidate plate nodes (j0 and its eight neighbours, live and of j0's crust type; n x 9, invalid entries j0) and
        their normalised weights."""
        g = self.g
        if not hasattr(self, '_tent_h0'):
            ax = g.axis & (g.nbr >= 0)
            self._tent_h0 = np.where(ax, g.dist, np.inf).min(1)
        nb = g.nbr[j0]
        cand = np.concatenate([j0[:, None], nb], 1)
        ok = cand >= 0
        cs = np.where(ok, cand, j0[:, None])
        ok &= pl.mask[cs] & (pl.cont[cs] == pl.cont[j0][:, None])
        q = g.xyz[ws]@pl.Q                                # world positions in the plate's frame (as _claims)
        d = 2*np.arcsin(np.minimum(np.linalg.norm(g.xyz[cs]-q[:, None, :], axis=2)/2, 1.))*g.radius   # arc, as g.dist
        ell = np.concatenate([self._tent_h0[j0][:, None], np.where(nb >= 0, g.dist[j0], np.inf)], 1)
        wt = np.where(ok, np.clip(1.-d/ell, 0., None), 0.)
        wt[:, 0] = np.maximum(wt[:, 0], 1e-12)            # the nearest node always carries some weight
        wt /= wt.sum(1, keepdims=True)
        return np.where(ok, cs, j0[:, None]), wt


    def _overlap_weights(self, pl, ws, j0):
        """remap_rule 'overlap' (pass 108, see PARAMS): each world cell sampled at nine points (its node and the points a
        third of the way to its eight neighbours), each point read from the plate node nearest it in the plate's frame
        (live and of j0's crust type, else j0); weight 1/9 per point: the cell's overlap with the plate cells."""
        g = self.g
        if not hasattr(self, '_overlap_sub'):
            nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
            sub = np.concatenate([g.xyz[:, None, :], g.xyz[:, None, :]+(g.xyz[nb]-g.xyz[:, None, :])/3.], 1)
            self._overlap_sub = sub/np.linalg.norm(sub, axis=2, keepdims=True)      # N x 9 x 3
        q = self._overlap_sub[ws].reshape(-1, 3)@pl.Q
        cs = fields.nearest_nodes(g, q).reshape(ws.size, 9)
        ok = pl.mask[cs] & (pl.cont[cs] == pl.cont[j0][:, None])
        cs = np.where(ok, cs, j0[:, None])
        return cs, np.full(cs.shape, 1./9.)


    def _tent_gather(self, idx, slot, gap, wthick, wsed, wrock, wborn, wtherm):
        """remap_rule 'tent': the world fields interpolated from each owner's plate raster; the weights are kept for the
        return transfer (self._tent[k] = (world nodes, candidates, weights))."""
        out = [v.copy() for v in (wthick, wsed, wrock, wborn, wtherm)]
        self._tent = {}
        self._tent_type = {}
        for k, pl in enumerate(self.plates):
            ws = np.flatnonzero((slot == k) & ~gap)
            if not ws.size:
                self._tent[k] = (ws, np.zeros((0, 9), np.int64), np.zeros((0, 9)))
                self._tent_type[k] = np.zeros(0, bool)
                continue
            weights = self._overlap_weights if self.p['remap_rule'] == 'overlap' else self._tent_weights
            cs, wt = weights(pl, ws, idx[k, ws])
            self._tent[k] = (ws, cs, wt)
            self._tent_type[k] = pl.cont[idx[k, ws]].copy()   # the crust type each world node read (its candidates')
            for o, src in zip(out[:4], (pl.thick, pl.sed, pl.rock, pl.born)):
                o[ws] = (src[cs].astype(np.float64)*wt).sum(1)
            th = pl.therm[cs].astype(np.float64)
            fin = np.isfinite(th)
            wf = np.where(fin, wt, 0.)
            sw = wf.sum(1)
            use = np.isfinite(th[:, 0]) & (sw > 0)        # the nearest node's validity decides; finite neighbours only
            out[4][ws[use]] = (np.where(fin, th, 0.)*wf).sum(1)[use]/sw[use]
        return tuple(out)


    def _tent_return(self, k, pl, wn, own, d, name):
        """remap_rule 'tent' / 'overlap': the per-area increment d (world) as per-area increments on plate pl's nodes, by
        the composition's weights (pass 112, journal 153: the material-aware return).

        A cell's increment goes to the plate nodes it read (its contributors), owned or not, as under 'nearest' (where
        the read nodes take the increments): a node the world shows is the material the world changed. Each increment is
        split into its gain and its loss before it is mixed, so a cut next to a fill cannot cancel and then be amplified (pass 109). Gains go
        back by the weights (each plate node the weighted mean of the gains over it). A cell's loss is debited from its
        contributors in proportion to weight x the material each has above its floor (`_return_floor`: oceanic
        thickness keeps its basement, `ocean_thick_m`, the floor the forearc cut respects; everything else 0), at most
        all of that material, so a contributor with nothing above its floor takes nothing (pass 111: the cut on a cell
        mixing an arc and basement went onto the basement). A live owned node no weight reaches takes the world node
        over it where the crust type agrees, its loss cut to its own material. Then each plate's gains and its losses,
        per crust type (world nodes typed by the crust they read), are closed to what its nodes can take by
        rescaling only the nodes that took that sign: the factor corrects the overlap areas' mismatch and stays near one,
        and no uniform correction moves untouched nodes (pass 109's did: journal 151 sec. 3). Booked, never spread: a
        cell's gain with no contributor (erosion_totals['return_unplaced_<name>_m3'], m^3), and a loss beyond its
        contributors' material or beyond a node's after the rescaling (['return_shortfall_<name>_m3'], loss not
        removed). (Pass 112's first version gave increments only to owned contributors and closed the whole world
        volume: a cell's gain on a hidden node, divided by the denormal gains of the two owned nodes of its type, made the
        factor infinite, Transfer Sweep seed 1004 at -134 Myr; renormalised over owned contributors, 90% of the loss it
        then booked was on cells whose contributors were all hidden or over new crust: journal 153.) Returns (plate
        nodes, increments)."""
        g = self.g
        a = g.area
        ws, cs, f = self._tent[k]
        wtype = self._tent_type[k]
        live = pl.mask & own & ~self.new_crust[wn]       # owned: the nodes no weight reaches take the world node over them
        take = pl.mask                                    # contributors: every plate node a world cell read
        fl = f*take[cs]
        sf = fl.sum(1)
        has = sf > 0
        fl = fl/np.where(has, sf, 1.)[:, None]            # weights over the contributors of each cell
        wt = fl*a[ws][:, None] if self.p['remap_rule'] == 'overlap' else fl   # overlap: the cell's area times its share
        avail = np.maximum(self._return_cur(pl, name).astype(np.float64)-self._return_floor(pl, name), 0.)
        dw = d[ws]
        gain, loss = np.maximum(dw, 0.), np.minimum(dw, 0.)
        ebar = (fl*avail[cs]).sum(1)                      # each cell's material above the floor
        ok = has & (ebar > 0)
        lw = np.where(ok, np.maximum(loss/np.where(ok, ebar, 1.), -1.), 0.)   # loss per unit of material, >= -1
        Lg = np.where(has, gain, 0.)                      # what the contributors can take
        Ll = lw*ebar
        W = np.bincount(cs.ravel(), weights=wt.ravel(), minlength=g.size)
        Sp = np.bincount(cs.ravel(), weights=(wt*gain[:, None]).ravel(), minlength=g.size)
        Sm = np.bincount(cs.ravel(), weights=(wt*lw[:, None]).ravel(), minlength=g.size)
        reach = take & (W > 0)
        cw = self.world['cont']
        rest = live & (W <= 0) & (cw[wn] == pl.cont)      # the world node over it, where the crust type agrees
        J = np.flatnonzero(reach | rest)
        r = reach[J]
        WJ = np.where(r, W[J], 1.)
        eJ = avail[J]
        dp = np.where(r, Sp[J]/WJ, np.maximum(d[wn[J]], 0.))
        dm = np.where(r, eJ*Sm[J]/WJ, np.maximum(np.minimum(d[wn[J]], 0.), -eJ))
        aJ = a[J]
        aw = a[ws]
        moved = short = unplaced = 0.
        for typ in (True, False):                         # each crust type, gains and losses each closed on their own
            mj = pl.cont[J] == typ
            mw = wtype == typ
            Vp, Vm = float((aw[mw]*gain[mw]).sum()), float((aw[mw]*loss[mw]).sum())
            Lp, Lm = float((aw[mw]*Lg[mw]).sum()), float((aw[mw]*Ll[mw]).sum())
            Pp, Pm = float((aJ[mj]*dp[mj]).sum()), float((aJ[mj]*dm[mj]).sum())
            fp = Lp/Pp if Pp > 0 else np.inf
            if np.isfinite(fp):
                moved += abs(Lp-Pp)
                dp[mj] *= fp
                unplaced += Vp-Lp
            else:
                dp[mj] = 0.
                unplaced += Vp
            fm = Lm/Pm if Pm < 0 else np.inf
            if np.isfinite(fm):
                moved += abs(Lm-Pm)
                want = dm[mj]*fm
                got = np.maximum(want, -eJ[mj])           # never more than the node's material
                short += float((aJ[mj]*(got-want)).sum())+(Lm-Vm)
                dm[mj] = got
            else:
                dm[mj] = 0.
                short -= Vm
        ET = self.erosion_totals
        for key, v in ((f'return_rescaled_{name}_m3', moved), (f'return_shortfall_{name}_m3', short),
                       (f'return_unplaced_{name}_m3', unplaced)):
            ET[key] = ET.get(key, 0.)+v
        return J, dp+dm


    def _return_cur(self, pl, name):
        return {'thick': pl.thick, 'sed': pl.sed, 'rock': pl.rock}[name]


    def _return_floor(self, pl, name):
        """The level a return's loss may not take a plate node below: oceanic crust keeps its basement (the forearc cut
        and the shelf never remove it), all else may go to zero."""
        if name == 'thick':
            return np.where(pl.cont, 0., float(self.p['ocean_thick_m']))
        return 0.


    def _match_reads(self, idx, slot, gap):
        """compose_map 'matched' (pass 107, see PARAMS): in place on idx, after ownership. Where world nodes of one
        plate read the same plate node, the nearest keeps it and each other reads the nearest live node of that plate
        within two rings of the shared node that no world node of the plate reads (four passes; a world node with no
        such node keeps the shared one)."""
        g = self.g
        for k, pl in enumerate(self.plates):
            ws = np.flatnonzero((slot == k) & ~gap)
            if ws.size == 0:
                continue
            loc = g.xyz[ws]@pl.Q                          # world nodes in the plate's frame
            for _ in range(4):
                js = idx[k, ws]
                d = ((loc-g.xyz[js])**2).sum(1)
                order = np.lexsort((d, js))
                first = np.ones(order.size, bool)
                first[1:] = js[order][1:] != js[order][:-1]
                dsp = order[~first]                       # positions in ws of displaced readers
                if dsp.size == 0:
                    break
                read = np.zeros(g.size, bool); read[js] = True
                n1 = g.nbr[js[dsp]]                       # the shared node's two rings
                n2 = g.nbr[np.where(n1 >= 0, n1, js[dsp][:, None])].reshape(dsp.size, -1)
                cand = np.concatenate([n1, n2], 1)
                ok = cand >= 0
                cs = np.where(ok, cand, 0)
                free = ok & pl.mask[cs] & ~read[cs]
                dc = np.where(free, ((g.xyz[cs]-loc[dsp][:, None, :])**2).sum(-1), np.inf)
                b = np.argmin(dc, 1)
                r = np.arange(dsp.size)
                has = np.isfinite(dc[r, b])
                if not has.any():
                    break
                dsp, ch, dd = dsp[has], cs[r, b][has], dc[r, b][has]
                o2 = np.lexsort((dd, ch))                 # one reader per freed node: the nearest
                f2 = np.ones(o2.size, bool)
                f2[1:] = ch[o2][1:] != ch[o2][:-1]
                take = o2[f2]
                idx[k, ws[dsp[take]]] = ch[take]


    def _spread_collision(self, dVc):
        """collision_spread 'earth' (see PARAMS): each receiving world node's collision volume (m^3) goes to the
        continental nodes of its own plate whose nearest receiving node it is, weighted by area x Earth's thickening
        kernel at their distance from it (graph distance through that plate's continental crust). Volume conserved;
        a receiving node with no such node keeps its volume."""
        g, w = self.g, self.world
        src = dVc > 0
        recv = w['cont'] | src
        ok = g.nbr >= 0
        nb = np.where(ok, g.nbr, 0)
        same = ok & recv[:, None] & recv[nb] & (self.owner[:, None] == self.owner[nb])
        cost = np.where(same, g.dist, 1e30)
        lab, d = nearest_label(np.where(src, np.arange(g.size), -1).astype(np.int64), g.nbr, cost)
        wt = np.where(recv & (lab >= 0) & (d < 1e29), g.area*collision_kernel(d), 0.)
        use = wt > 0
        W = np.bincount(lab[use], weights=wt[use], minlength=g.size)
        out = np.zeros(g.size)
        has = src & (W > 0)
        sel = use & has[np.where(lab >= 0, lab, 0)]
        np.add.at(out, np.flatnonzero(sel), dVc[lab[sel]]*wt[sel]/W[lab[sel]])
        out[src & ~has] += dVc[src & ~has]
        self.erosion_totals['collision_spread_m3'] = self.erosion_totals.get('collision_spread_m3', 0.)+float(
            np.abs(out-dVc).sum()/2)
        return out


    def _return_pairs(self, pl, wn, own, ws, js, nread):
        """return_transfer 'conservative' / 'matched' (pass 107, see PARAMS): the plate nodes J that take a world
        increment this step and the world node wJ each takes it from (each plate node once).

        'conservative': every live plate node the plate still owns (not this step's new crust), from the world node over
        it (wn, the forward map). 'matched': each plate node a world node reads keeps the increment of its reader (the
        nearest reader to it when several read it), so the world node finds its own change on the material it read; the
        live owned nodes no world node reads take the increment of a reader displaced that way (the nearest among the
        world node over them and its eight neighbours), or else of the world node over them."""
        g = self.g
        live = np.flatnonzero(pl.mask & own & ~self.new_crust[wn])
        if self.p['return_transfer'] == 'conservative':
            return live, wn[live]
        if ws.size == 0:
            return live, wn[live]
        loc = g.xyz@pl.Q.T                               # plate nodes' world positions
        d = np.einsum('ij,ij->i', g.xyz[ws]-loc[js], g.xyz[ws]-loc[js])
        order = np.lexsort((d, js))                      # by plate node, nearest reader first
        first = np.ones(order.size, bool)
        first[1:] = js[order][1:] != js[order][:-1]
        prim = order[first]                              # one reader per read plate node
        disp = np.zeros(g.size, bool)
        disp[ws[order[~first]]] = True                   # readers displaced from a shared plate node
        un = live[nread[live] == 0]
        src = wn[un]
        if un.size and disp.any():
            cand = np.concatenate([src[:, None], g.nbr[src]], 1)
            ok = cand >= 0
            cs = np.where(ok, cand, 0)
            dc = ((g.xyz[cs]-loc[un][:, None, :])**2).sum(-1)
            dc = np.where(ok & disp[cs], dc, np.inf)
            b = np.argmin(dc, 1)
            has = np.isfinite(dc[np.arange(un.size), b])
            src = np.where(has, cs[np.arange(un.size), b], src)
        return np.concatenate([js[prim], un]), np.concatenate([ws[prim], src])


    def _conserved(self, dj, aj, dw, aw, name):
        """return_transfer 'conservative': per-area increments dj on plate nodes of area aj (read from the world node
        over each) rescaled so that their gains and their losses each carry the world volume of the nodes the plate
        owns (dw on areas aw). A part the plate nodes do not see at all (no gain where the world has one) is spread
        uniformly by area. The rescaling's volume is recorded in erosion_totals['return_rescaled_<name>_m3']."""
        dj = np.asarray(dj, np.float64)
        dw = np.asarray(dw, np.float64)
        out = dj.copy()
        rest = 0.
        for sg in (1., -1.):
            mj = sg*dj > 0
            Pj = float((aj*dj)[mj].sum())
            Pw = float((aw*dw)[sg*dw > 0].sum())
            if Pj != 0.:
                out[mj] = dj[mj]*(Pw/Pj)
            else:
                rest += Pw
        if rest != 0.:
            out += rest/aj.sum()
        ET = self.erosion_totals
        key = f'return_rescaled_{name}_m3'
        ET[key] = ET.get(key, 0.)+float((np.abs(out-dj)*aj).sum())
        return out


    def _forward_nodes(self, pl):
        """World node nearest to each of plate pl's nodes (the forward map, plate node j -> world node nearest Q j). A
        pure function of the rotation, so it is kept for the step by the rotation's bytes (journal 154: _opened_area,
        _update_rasters and a split read the same map); the array is read-only."""
        cache = self.__dict__.setdefault('_forward_cache', {})
        key = pl.Q.tobytes()
        wn = cache.get(key)
        if wn is None:
            wn = fields.nearest_nodes(self.g, self.g.xyz@pl.Q.T)
            wn.flags.writeable = False
            cache[key] = wn
        return wn


    def _update_rasters(self, dthick, idx, dsed=None, drock=None):
        """Update every plate raster after a step.

        Ownership: forward lookup (plate node j -> world node). Overlapped crust stays hidden for hidden_steps before
        it is consumed (subducted ocean, or continental crust whose volume is added to the overriding continent:
        collision). Thickness changes go to the plate nodes the world actually reads (backward samples), so every
        change feeds back: diffusion increments (per area) are averaged over the world nodes reading a plate node;
        collision volume is added as volume. Plate nodes the world does not read keep their values.

        return_transfer 'conservative' (pass 107, see PARAMS): the per-area increments (thickness, sediment, rock),
        orogeny marks, thermal-age changes and arcs turned continental go instead to every live plate node the plate
        still owns, from the world node over it (the forward map), with each plate's gains and losses rescaled to the
        world volume of the nodes it owns; collision volume and new crust as above."""
        g, w, p = self.g, self.world, self.p
        if p['return_transfer'] not in ('read', 'conservative', 'matched', 'state'):
            raise ValueError(f"return_transfer {p['return_transfer']!r}: 'read', 'conservative', 'matched' or 'state'")
        tent = p['remap_rule'] in ('tent', 'overlap')
        if tent and p['return_transfer'] != 'read':
            raise ValueError(f"remap_rule {p['remap_rule']!r} needs return_transfer 'read'")
        a = g.area
        wns, owns, consumed, booked = [], [], [], []
        dVc = np.zeros(g.size)                          # collision volume arriving at world nodes (m^3)
        collided = np.zeros(g.size, bool)
        subducted = np.zeros(g.size, bool)
        for k, pl in enumerate(self.plates):            # pass 1: hidden / consumed crust
            wn = self._forward_nodes(pl)
            own = self.owner_slot[wn] == k
            hidden = pl.mask & ~own
            pl.under[own] = 0
            pl.under[hidden] += 1
            gone = hidden & (pl.under > p['hidden_steps'])
            cc = gone & pl.cont & w['cont'][wn]
            np.add.at(dVc, wn[cc], pl.thick[cc].astype(np.float64)*a[cc])
            collided[wn[cc]] = True
            SL = self.sed_ledger
            if p['continental_sediment']:                # its sediment is folded into the receiving belt's crust
                np.add.at(dVc, wn[cc], pl.sed[cc].astype(np.float64)*a[cc])
                self.erosion_totals['collision_sediment_to_crust_m3'] = self.erosion_totals.get(
                    'collision_sediment_to_crust_m3', 0.)+float((pl.sed[cc]*a[cc]).sum())
            SL['collision_out'] -= float((pl.sed[cc]*a[cc]).sum())
            SL['cont_subducted'] -= float((pl.sed[gone & pl.cont & ~w['cont'][wn]]*a[gone & pl.cont & ~w['cont'][wn]]).sum())
            subducted[wn[gone & ~pl.cont]] = True
            if p['trench_flux'] == 'consumed':           # B-d2: the oceanic area consumed feeds next step's trench fluxes
                if getattr(self, 'consumed_pending', None) is None:
                    self.consumed_pending = np.zeros(g.size)
                oc_all = gone & ~pl.cont
                np.add.at(self.consumed_pending, wn[oc_all], a[oc_all])
            if p['arc_accretion']:                        # B-d: buoyant arc crust docks on the overriding plate
                oc = gone & ~pl.cont & (pl.thick > p['ocean_thick_m'])
                ex = pl.thick[oc].astype(np.float64)-p['ocean_thick_m']
                acc = np.where(pl.thick[oc] >= p['arc_buoyant_thick_m'], ex, np.minimum(ex, p['arc_accrete_upper_m']))
                on_c = w['cont'][wn[oc]]
                np.add.at(dVc, wn[oc][on_c], (acc*a[oc])[on_c])  # on continental crust: collision volume, as now
                # on an oceanic overriding plate the docked crust joins that plate's arc: it is spread next step over
                # the arc band (arc_d0_km +- arc_sigma_km) of the nearest overriding trench node, volume conserved
                # (placed on the single receiving node it made one-node spikes that turned continental as specks)
                if not hasattr(self, 'dock_pending'):
                    self.dock_pending = np.zeros(g.size)
                np.add.at(self.dock_pending, wn[oc][~on_c], (acc*a[oc])[~on_c])
                collided[wn[oc][on_c]] = True
                ET = self.erosion_totals
                ET['arc_accreted_m3'] = ET.get('arc_accreted_m3', 0.)+float((acc*a[oc]).sum())
                ET['arc_accreted_on_continent_m3'] = ET.get('arc_accreted_on_continent_m3', 0.)+float((acc*a[oc])[on_c].sum())
            sg = gone & ~pl.cont & (pl.sed > 0)          # trench sediment: part is scraped onto the overriding plate
            np.add.at(dVc, wn[sg], p['sediment_accretion']*pl.sed[sg].astype(np.float64)*a[sg])
            self.erosion_totals['sediment_accreted_m3'] = self.erosion_totals.get('sediment_accreted_m3', 0.)+float(
                (p['sediment_accretion']*pl.sed[sg]*a[sg]).sum())
            self.erosion_totals['sediment_subducted_m3'] = self.erosion_totals.get('sediment_subducted_m3', 0.)+float(
                ((1-p['sediment_accretion'])*pl.sed[sg]*a[sg]).sum())
            L = self.ledger
            L['collision_out'] -= float((pl.thick[cc]*a[cc]).sum())
            lost = gone & pl.cont & ~w['cont'][wn]
            self.area_ledger['collision'] -= float(a[cc].sum())
            self.area_ledger['subduction'] -= float(a[lost].sum())
            L['cont_subducted'] -= float((pl.thick[lost]*a[lost]).sum())
            wns.append(wn); owns.append(own); consumed.append(gone)
            g1 = gone & pl.cont                          # booked above at their values now (pass 2 may still change them)
            booked.append((cc, np.where(g1, pl.thick, 0.).astype(np.float64), np.where(g1, pl.sed, 0.).astype(np.float64)))
        self.collided, self.subducted = collided, subducted
        if p['collision_spread'] == 'earth' and dVc.any():
            dVc = self._spread_collision(dVc)
        elif p['collision_spread'] not in ('node', 'earth'):
            raise ValueError(f"collision_spread {p['collision_spread']!r}: 'node' or 'earth'")
        if hasattr(self, 'tectonic'):
            self.tectonic = self.tectonic+dVc/a          # collision volume, as thickness at the receiving world node
        marked = collided | self._orogen_marks()
        self.orogen_marked = marked
        w['orog'] = np.where(marked, self.t, w['orog'])
        thick_l = np.zeros(g.size, bool)
        if p['collision_lithosphere'] == 'thicken' and p['lithosphere'] == 'thermal_age':
            # C, pure shear (journal 101): the mantle lithosphere thickens with the crust that collision thickens,
            # L -> L (T + dT) / T, the mirror of the stretching rule; no reset at collision contacts
            w['therm'], thick_l = thickened_therm(w['therm'], w['thick'], np.where(w['cont'], dVc/a, 0.), self.t, p)
        elif p['collision_lithosphere'] not in ('reset', 'thicken'):
            raise ValueError(f"collision_lithosphere {p['collision_lithosphere']!r}: 'reset' or 'thicken'")
        reset = self._thermal_resets(collided)           # C: collision contacts ('reset') and continental arcs
        w['therm'][reset] = self.t
        tchg = self.therm_changed | reset | thick_l
        for k, pl in enumerate(self.plates):            # pass 2: ownership, new crust, thickness changes
            wn, own, gone = wns[k], owns[k], consumed[k]
            # a plate node switches on for genuinely new crust (a ridge gap), or for oceanic crust the plate already
            # covers through a neighbouring node (filling the staircase behind a spreading ridge). Continental crust
            # is never duplicated that way (the forward and backward nearest-node maps are not exact inverses), so
            # the continental inventory is conserved.
            new = own & ~pl.mask & (self.new_crust[wn] | ~w['cont'][wn])
            ws = np.flatnonzero((self.owner_slot == k) & ~self.new_crust)   # world nodes read from this plate
            js = idx[k, ws]
            nread = np.bincount(js, minlength=g.size)
            read = nread > 0
            before_v = float((pl.thick*a*(pl.mask & pl.cont)).sum())
            vs = self._plate_cont_sed(pl)
            # return_transfer 'conservative' / 'matched': plate nodes J and the world node wJ whose increment each takes
            J, wJ = (self._return_pairs(pl, wn, own, ws, js, nread) if p['return_transfer'] in ('conservative', 'matched')
                     else (None, None))
            cons = J is not None and J.size > 0
            if cons:
                d_vol = np.bincount(js, weights=dVc[ws], minlength=g.size)/a    # collision volume as 'read'
                pl.thick[read] += d_vol[read].astype(np.float32)
                pl.thick[J] += self._conserved(dthick[wJ], a[J], dthick[ws], a[ws], 'thick').astype(np.float32)
                if dsed is not None:
                    s1 = pl.sed[J]+self._conserved(dsed[wJ], a[J], dsed[ws], a[ws], 'sed')
                    self.erosion_totals['return_sed_clamped_m3'] = self.erosion_totals.get('return_sed_clamped_m3', 0.)+float(
                        (np.maximum(-s1, 0.)*a[J]).sum())
                    pl.sed[J] = np.maximum(s1, 0).astype(np.float32)
                if drock is not None:
                    pl.rock[J] = np.clip(pl.rock[J]+self._conserved(drock[wJ], a[J], drock[ws], a[ws], 'rock'), 0.,
                                         pl.sed[J]).astype(np.float32)
            elif tent:                                  # remap_rule 'tent': the composition's weights, conserved
                d_vol = np.bincount(js, weights=dVc[ws], minlength=g.size)/a    # collision volume as 'read'
                pl.thick[read] += d_vol[read].astype(np.float32)
                Jt, dt_ = self._tent_return(k, pl, wn, own, dthick, 'thick')
                t1 = pl.thick[Jt]+dt_                     # a node the area-mean view hides can be thinned below zero
                self.erosion_totals['return_thick_clamped_m3'] = self.erosion_totals.get('return_thick_clamped_m3', 0.)+float(
                    (np.maximum(-t1, 0.)*a[Jt]).sum())     # (pass 109): clamped at zero like the sediment, booked
                pl.thick[Jt] = np.maximum(t1, 0.).astype(np.float32)
                if dsed is not None:
                    Js, ds_ = self._tent_return(k, pl, wn, own, dsed, 'sed')
                    s1 = pl.sed[Js]+ds_
                    self.erosion_totals['return_sed_clamped_m3'] = self.erosion_totals.get('return_sed_clamped_m3', 0.)+float(
                        (np.maximum(-s1, 0.)*a[Js]).sum())
                    pl.sed[Js] = np.maximum(s1, 0).astype(np.float32)
                if drock is not None:
                    Jr, dr_ = self._tent_return(k, pl, wn, own, drock, 'rock')
                    pl.rock[Jr] = np.clip(pl.rock[Jr]+dr_, 0., pl.sed[Jr]).astype(np.float32)
            else:
                d_int = np.bincount(js, weights=dthick[ws], minlength=g.size)/np.maximum(nread, 1)
                d_vol = np.bincount(js, weights=dVc[ws], minlength=g.size)/a
                pl.thick[read] += (d_int[read]+d_vol[read]).astype(np.float32)
                if dsed is not None:
                    d_sed = np.bincount(js, weights=dsed[ws], minlength=g.size)/np.maximum(nread, 1)
                    pl.sed[read] = np.maximum(pl.sed[read]+d_sed[read], 0).astype(np.float32)
                if drock is not None:                       # sediment_rock: the rock layer's change, as the sediment's
                    d_rk = np.bincount(js, weights=drock[ws], minlength=g.size)/np.maximum(nread, 1)
                    pl.rock[read] = np.clip(pl.rock[read]+d_rk[read], 0., pl.sed[read]).astype(np.float32)
            self.ledger['increments'] += float((pl.thick*a*(pl.mask & pl.cont)).sum())-before_v
            vs1 = self._plate_cont_sed(pl)
            self.sed_ledger['increments'] += vs1-vs
            if p['return_transfer'] == 'state':
                # 'state': the live owned nodes no world node read take the world's end state at the node over them
                # (the surface the world showed there), where the crust type agrees; booked as 'return_state'
                un = np.flatnonzero(pl.mask & own & ~self.new_crust[wn] & (nread == 0))
                un = un[pl.cont[un] == w['cont'][wn[un]]]
                if un.size:
                    wu = wn[un]
                    v0 = float((pl.thick*a*(pl.mask & pl.cont)).sum())
                    pl.thick[un] = w['thick'][wu]
                    pl.sed[un] = np.maximum(w['sed'][wu], 0.)
                    pl.rock[un] = np.clip(w['rock'][wu], 0., pl.sed[un])
                    pl.therm[un] = w['therm'][wu]
                    pl.orog[un] = w['orog'][wu]
                    self.ledger['return_state'] = self.ledger.get('return_state', 0.)+float(
                        (pl.thick*a*(pl.mask & pl.cont)).sum())-v0
                    vs2 = self._plate_cont_sed(pl)
                    self.sed_ledger['return_state'] = self.sed_ledger.get('return_state', 0.)+vs2-vs1
                    vs1 = vs2
            nc = getattr(self, 'new_continent', None)
            if cons:                                    # 'conservative': the plate nodes under the arcs (forward map)
                conv = nc[wJ] if nc is not None else np.zeros(J.size, bool)
                cj = J
            else:
                conv = nc[ws] if nc is not None else np.zeros(ws.size, bool)
                cj = js
            if conv.any():                              # each plate node once: two world nodes can read the same one
                uc = np.unique(cj[conv])                   # (booked twice until the thirty-ninth pass; ledger only)
                pl.cont[uc] = True
                pl.orog[uc] = self.t
                self.ledger['arc_continent'] = self.ledger.get('arc_continent', 0.)+float((pl.thick[uc]*a[uc]).sum())
                self.sed_ledger['arc_continent'] += self._plate_cont_sed(pl)-vs1
            hit = np.zeros(g.size, bool)
            if cons:                                    # 'conservative': marks and thermal ages by the forward map
                hit[J[marked[wJ]]] = True
                pl.orog[hit] = self.t
                sel = tchg[wJ]
                pl.therm[J[sel]] = w['therm'][wJ[sel]]
            else:
                hit[js[marked[ws]]] = True
                pl.orog[hit] = self.t
                sel = tchg[ws]                           # C: thermal ages changed at world nodes this plate feeds
                pl.therm[js[sel]] = w['therm'][ws[sel]]
            pl.cont[new] = w['cont'][wn[new]]
            pl.thick[new] = w['thick'][wn[new]]
            pl.born[new] = w['born'][wn[new]]
            pl.orog[new] = w['orog'][wn[new]]
            pl.therm[new] = w['therm'][wn[new]]
            pl.sed[new] = w['sed'][wn[new]]
            pl.rock[new] = w['rock'][wn[new]]
            pl.under[new] = 0
            pl.slab[new] = 0
            self.ledger['new_nodes'] += float((pl.thick[new]*a[new]*pl.cont[new]).sum())
            self.sed_ledger['new_nodes'] += float((pl.sed[new]*a[new]*pl.cont[new]).sum())
            # consumed nodes can still be read by a world node of this plate (the forward and backward maps are not
            # exact inverses) and get this step's increments, or become continental (arcs), after pass 1 booked them:
            # book what leaves with them at its value now (thirty-eighth pass; ledgers only, the fields are unchanged)
            cc1, th1, sd1 = booked[k]
            now_c = gone & pl.cont
            d_th = (np.where(now_c, pl.thick, 0.).astype(np.float64)-th1)*a
            d_sd = (np.where(now_c, pl.sed, 0.).astype(np.float64)-sd1)*a
            self.ledger['collision_out'] -= float(d_th[cc1].sum())
            self.ledger['cont_subducted'] -= float(d_th[~cc1].sum())
            self.sed_ledger['collision_out'] -= float(d_sd[cc1].sum())
            self.sed_ledger['cont_subducted'] -= float(d_sd[~cc1].sum())
            pl.mask[gone] = False
            pl.mask[new] = True
            pl.cont[~pl.mask] = False
            pl.slab[~pl.mask] = 0
            hole, src = plate_holes(pl.mask, g.nbr)
            if hole.any():
                hs = src[hole]
                pl.mask[hole] = True
                for arr in (pl.cont, pl.thick, pl.born, pl.orog, pl.sed, pl.slab, pl.therm, pl.rock):
                    arr[hole] = arr[hs]
                self.ledger['holes'] += float((pl.thick[hole]*a[hole]*pl.cont[hole]).sum())
                self.area_ledger['holes'] = self.area_ledger.get('holes', 0.)+float(a[hole & pl.cont].sum())
                self.sed_ledger['holes'] += float((pl.sed[hole]*a[hole]*pl.cont[hole]).sum())
            self.holes_filled = getattr(self, 'holes_filled', 0)+int(hole.sum())
        # collision volume that reached a world node no plate node reads is lost. Its removal is already booked in
        # pass 1 (collision_out) and what arrives is booked as increments, so it is counted here, not booked again
        # (until the thirty-eighth pass it was also booked as ledger 'cont_lost_other': a double count)
        self.erosion_totals['collision_volume_unread_m3'] = self.erosion_totals.get('collision_volume_unread_m3', 0.)+float(
            dVc[self.new_crust].sum())
        self._wn = {pl.id: wn for pl, wn in zip(self.plates, wns)}   # raster -> world maps (B-b reuses them)
