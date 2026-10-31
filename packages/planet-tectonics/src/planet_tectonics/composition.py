"""PlateSim, composition: which plate owns each world node (claims, scores, ocean-ocean polarity), the world fields,
and the plate-raster update after a step (hidden and consumed crust, collision volume, increments, ledgers)."""
import numpy as np

from planet_grid import fields

from planet_tectonics.graphops import nearest_label, plate_holes


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
        wcont = np.where(gap, False, cont_k[slot, cols])
        wthick = np.where(gap, p['ocean_thick_m'], np.stack([pl.thick[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        wborn = np.where(gap, self.t, np.stack([pl.born[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        worog = np.where(gap, np.nan, np.stack([pl.orog[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        wsed = np.where(gap, 0., np.stack([pl.sed[idx[k]] for k, pl in enumerate(self.plates)])[slot, cols]).astype(np.float64)
        self.new_crust = gap
        self.world = {'cont': wcont, 'thick': wthick, 'born': wborn, 'orog': worog, 'sed': wsed}


    def _update_rasters(self, dthick, idx, dsed=None):
        """Update every plate raster after a step.

        Ownership: forward lookup (plate node j -> world node). Overlapped crust stays hidden for hidden_steps before
        it is consumed (subducted ocean, or continental crust whose volume is added to the overriding continent:
        collision). Thickness changes go to the plate nodes the world actually reads (backward samples), so every
        change feeds back: diffusion increments (per area) are averaged over the world nodes reading a plate node;
        collision volume is added as volume. Plate nodes the world does not read keep their values."""
        g, w, p = self.g, self.world, self.p
        a = g.area
        wns, owns, consumed, booked = [], [], [], []
        dVc = np.zeros(g.size)                          # collision volume arriving at world nodes (m^3)
        collided = np.zeros(g.size, bool)
        subducted = np.zeros(g.size, bool)
        for k, pl in enumerate(self.plates):            # pass 1: hidden / consumed crust
            wn = fields.nearest_nodes(g, g.xyz@pl.Q.T)
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
        if hasattr(self, 'tectonic'):
            self.tectonic = self.tectonic+dVc/a          # collision volume, as thickness at the receiving world node
        marked = collided | self._orogen_marks()
        self.orogen_marked = marked
        w['orog'] = np.where(marked, self.t, w['orog'])
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
            d_int = np.bincount(js, weights=dthick[ws], minlength=g.size)/np.maximum(nread, 1)
            d_vol = np.bincount(js, weights=dVc[ws], minlength=g.size)/a
            pl.thick[read] += (d_int[read]+d_vol[read]).astype(np.float32)
            if dsed is not None:
                d_sed = np.bincount(js, weights=dsed[ws], minlength=g.size)/np.maximum(nread, 1)
                pl.sed[read] = np.maximum(pl.sed[read]+d_sed[read], 0).astype(np.float32)
            self.ledger['increments'] += float((pl.thick*a*(pl.mask & pl.cont)).sum())-before_v
            vs1 = self._plate_cont_sed(pl)
            self.sed_ledger['increments'] += vs1-vs
            conv = self.new_continent[ws] if hasattr(self, 'new_continent') else np.zeros(ws.size, bool)
            if conv.any():                              # each plate node once: two world nodes can read the same one
                uc = np.unique(js[conv])                   # (booked twice until the thirty-ninth pass; ledger only)
                pl.cont[uc] = True
                pl.orog[uc] = self.t
                self.ledger['arc_continent'] = self.ledger.get('arc_continent', 0.)+float((pl.thick[uc]*a[uc]).sum())
                self.sed_ledger['arc_continent'] += self._plate_cont_sed(pl)-vs1
            hit = np.zeros(g.size, bool); hit[js[marked[ws]]] = True
            pl.orog[hit] = self.t
            pl.cont[new] = w['cont'][wn[new]]
            pl.thick[new] = w['thick'][wn[new]]
            pl.born[new] = w['born'][wn[new]]
            pl.orog[new] = w['orog'][wn[new]]
            pl.sed[new] = w['sed'][wn[new]]
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
                for arr in (pl.cont, pl.thick, pl.born, pl.orog, pl.sed, pl.slab):
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
