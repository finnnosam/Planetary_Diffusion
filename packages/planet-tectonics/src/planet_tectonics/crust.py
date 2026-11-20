"""PlateSim, crust: rifting (node and pure-shear rules), sediment basins, subduction-zone thickening and trench fluxes,
orogeny marks and margin shortening."""
import numpy as np

from planet_landscape import graph_landscape as gl

from planet_tectonics.crust_kernels import margin_shortening, rift_pure_shear, spread_to_band, trench_flux_consumed
from planet_tectonics.graphops import nearest_label
from planet_tectonics.lithosphere import stretched_age
from planet_tectonics.params import BACKARC_PROFILE_KM, BACKARC_PROFILE_W, COLLISION, SUBDUCTION


class CrustMixin:
    """PlateSim, crust: rifting (node and pure-shear rules), sediment basins, subduction-zone thickening and trench fluxes,
    orogeny marks and margin shortening."""

    def _rift(self):
        """Gaps opening between continental crust become stretched continental crust (see PARAMS). Returns the
        thinning of the continental neighbours that supply the volume (m, <= 0)."""
        g, p, w = self.g, self.p, self.world
        out = np.zeros(g.size)
        gap = self.new_crust
        if not gap.any():
            return out
        rows = np.flatnonzero(gap)
        nb = g.nbr[rows]
        ok = nb >= 0
        nbs = np.where(ok, nb, 0)
        solid = ok & ~gap[nbs]
        cn = solid & w['cont'][nbs]
        ncont = cn.sum(1)
        if p['rift_rule'] == 'pure_shear':
            return self._rift_pure_shear(rows, ncont, solid, cn, nbs)
        meanT = np.where(cn, w['thick'][nbs], 0.).sum(1)/np.maximum(ncont, 1)
        Tg = p['rift_thinning']*meanT
        rift = (ncont > 0) & (ncont*2 >= solid.sum(1)) & (Tg >= p['breakup_thick_m'])
        if not rift.any():
            return out
        r, Tr, cr = rows[rift], Tg[rift], cn[rift]
        a = g.area
        if p['rift_width_km'] > 0:                       # necking zone: the volume comes from a tapered band
            src = np.zeros(g.size, bool); src[r] = True
            D = max(p['rift_width_km']*1e3, 1.5*float(g.dist[g.nbr >= 0].max()))   # at least the direct neighbours
            d = gl.distance_to(src, g.nbr, g.dist, D)
            band = w['cont'] & ~gap & (d < D)
            wt = np.where(band, (1-d/D)*a, 0.)
            take = (Tr*a[r]).sum()*wt/max(wt.sum(), 1e-9)
        else:
            share = (Tr*a[r]/ncont[rift])[:, None]*cr      # volume each continental neighbour gives (m^3)
            take = np.zeros(g.size)
            np.add.at(take, nbs[rift][cr], share[cr])
        out = -np.minimum(take/a, .5*w['thick'])
        w['cont'][r] = True
        w['thick'][r] = Tr
        w['born'][r] = np.where(cr, w['born'][nbs[rift]], np.inf).min(1)   # stretched old crust keeps its age
        w['therm'][r] = self.t                           # C: new rift crust, new lithosphere
        self.therm_changed[r] = True
        self.erosion_totals['rift_nodes'] = self.erosion_totals.get('rift_nodes', 0)+int(r.size)
        return out


    def _rift_pure_shear(self, rows, ncont, solid, cn, nbs):
        """rift_rule 'pure_shear' (see PARAMS): gaps between continental crust are the new area of a stretching zone."""
        g, p, w = self.g, self.p, self.world
        out = np.zeros(g.size)
        cand = (ncont > 0) & (ncont*2 >= solid.sum(1))
        if not cand.any():
            return out
        rc = rows[cand]
        zone = w['cont'] & ~self.new_crust
        if p['rift_width_km'] <= 0:
            raise ValueError("rift_rule 'pure_shear' needs rift_width_km > 0 (the stretching zone's half-width)")
        W = p['rift_width_km']*1e3
        dT = np.zeros(g.size)
        th = w['thick'].astype(np.float64)
        _, Tg, dT_s = rift_pure_shear(rc.astype(np.int64), g.area, th+dT, zone, g.nbr, g.dist, W, p['breakup_thick_m'])
        dT = dT+dT_s
        rift = Tg > 0.
        if not rift.any():
            return dT
        self.area_ledger['rift'] += float(g.area[rc[rift]].sum())
        r, cr = rc[rift], cn[cand][rift]
        w['cont'][r] = True
        w['thick'][r] = Tg[rift]
        w['born'][r] = np.where(cr, w['born'][nbs[cand][rift]], np.inf).min(1)   # stretched old crust keeps its age
        w['therm'][r] = self.t                           # C: new rift crust, new lithosphere
        self.therm_changed[r] = True
        self.erosion_totals['rift_nodes'] = self.erosion_totals.get('rift_nodes', 0)+int(r.size)
        self.erosion_totals['rift_stretched_m3'] = self.erosion_totals.get('rift_stretched_m3', 0.)+float((Tg[rift]*g.area[r]).sum())
        return dT


    def _sediment_basins_to_continent(self):
        """Oceanic nodes carrying more than basin_to_continent_sed_m of sediment, or sediment above sea level, become
        continental crust. Its thickness is set so the node keeps its present elevation (the isostatic equivalent)."""
        g, p, w = self.g, self.p, self.world
        dT = np.zeros(g.size); dS = np.zeros(g.size)
        if not p['basin_to_continent']:
            return dT, dS
        z = self.elevation()
        conv = (~w['cont']) & ((w['sed'] > p['basin_to_continent_sed_m']) | ((w['sed'] > 500.) & (z > 0)))
        if not conv.any():
            return dT, dS
        rm, rc, rw = p['rho_mantle'], p['rho_crust'], p['rho_water']
        zr = z-self.dyn
        dry = np.where(zr > 0, zr, zr*(rm-rw)/rm)
        T_new = p['cont_ref_thick_m']+(dry-p['cont_ref_elev_m'])*rm/(rm-rc)
        T_new = np.maximum(T_new, w['thick'])
        dT[conv] = T_new[conv]-w['thick'][conv]
        dS[conv] = -w['sed'][conv]
        a = g.area
        self.erosion_totals['basin_to_continent_m3'] = self.erosion_totals.get('basin_to_continent_m3', 0.)+float((w['sed'][conv]*a[conv]).sum())
        w['cont'] = w['cont'] | conv
        w['therm'][conv] = w['born'][conv]               # C: the floor's lithosphere keeps the floor's age
        self.therm_changed |= conv
        w['thick'] = w['thick']+dT
        w['sed'] = w['sed']+dS
        self.new_continent = self.new_continent | conv
        return dT, dS


    def _subduction_thickening(self):
        """Thickening (m over this step) of the overriding plate near trenches, from the previous step's boundaries."""
        g, p, w = self.g, self.p, self.world
        out = np.zeros(g.size)
        self.belt = np.zeros(g.size, bool)
        if p['belt_lithosphere'] == 'shorten':          # C: crust margin shortening adds to belt nodes (journal 105)
            self.belt_gain = np.zeros(g.size)
        if not hasattr(self, 'btype'):
            return out
        src = (self.btype == SUBDUCTION) & ~self.trench          # overriding side of a subduction boundary
        if not src.any():
            return out
        lab = np.where(src, np.arange(g.size), -1).astype(np.int64)
        near, d = nearest_label(lab, g.nbr, g.dist)
        ok = (near >= 0) & (d < 1000e3)
        ok[ok] &= self.owner[ok] == self.owner[near[ok]]           # stay on the overriding plate
        v = np.clip(self.conv_rate[np.clip(near, 0, None)], 0, None)  # m/yr closing
        dt = p['dt_myr']*1e6
        cont = w['cont']
        shorten = p['margin_shortening']
        H_cont = p['arc_H_cont_m'] if shorten else p['subduction_H_cont_m']
        terms = ((ok & cont, H_cont, p['subduction_d0_cont_km']*1e3, p['subduction_sigma_cont_km']*1e3),
                 (ok & ~cont, p['arc_H_ocean_m'], p['arc_d0_km']*1e3, p['arc_sigma_km']*1e3),
                 (ok, -p['subduction_erosion_H_m'], p['subduction_erosion_d0_km']*1e3, p['subduction_erosion_sigma_km']*1e3))
        mode = p['trench_flux']
        if mode == 'band':
            for mask, H, d0, sg in terms:
                if H == 0:
                    continue
                shape = np.exp(-.5*((d-d0)/sg)**2)/(sg*np.sqrt(2*np.pi))
                out = out+np.where(mask, v*H*shape*dt, 0.)
        elif mode == 'consumed':
            A = getattr(self, 'consumed_pending', None)
            self.consumed_pending = np.zeros(g.size)
            f, routed, unrouted = trench_flux_consumed(A, near, d, ok, terms, g.area, .5*np.sqrt(g.area.mean()),
                                                       self.owner, g.nbr, g.dist)
            out = out+f
            ET = self.erosion_totals
            ET['trench_flux_routed_m2'] = ET.get('trench_flux_routed_m2', 0.)+routed
            ET['trench_flux_unrouted_m2'] = ET.get('trench_flux_unrouted_m2', 0.)+unrouted
        else:
            raise ValueError(f"trench_flux {mode!r}: 'band' or 'consumed'")
        # erosion takes at most 30% of a continental column per step; on oceanic crust it removes only arc crust built
        # above the normal 7 km (an intra-oceanic forearc), never the oceanic basement: unbounded, it dug pits to
        # -10 km that trailed behind moving trenches (step 7 run)
        floor = np.where(cont, -.3*w['thick'], -np.clip(w['thick']-p['ocean_thick_m'], 0, None))
        out = np.maximum(out, floor)
        # the subduction belt on continental crust (orogen_rule 'belt'): within 2 sigma of d0, trench converging
        d0c, sgc = p['subduction_d0_cont_km']*1e3, p['subduction_sigma_cont_km']*1e3
        self.belt = ok & cont & (np.abs(d-d0c) <= 2*sgc) & (v > 0)
        for key, m in (('subduction_on_continents_m3', cont), ('subduction_on_ocean_arcs_m3', ~cont)):
            self.erosion_totals[key] = self.erosion_totals.get(key, 0.)+float((out*g.area)[m].sum())
        self.erosion_totals['subduction_added_m3'] = self.erosion_totals.get('subduction_added_m3', 0.)+float((out*g.area).sum())
        dock = getattr(self, 'dock_pending', None)   # B-d: arc crust docked on oceanic overriding plates last step
        if dock is not None and dock.any():          # (moved crust, not new: kept out of the arc-growth totals)
            dt_dock = spread_to_band(dock, near, d, ok & ~cont, p['arc_d0_km']*1e3, p['arc_sigma_km']*1e3, g.area)
            self.erosion_totals['arc_docked_spread_m3'] = self.erosion_totals.get('arc_docked_spread_m3', 0.)+float((dt_dock*g.area).sum())
            out = out+dt_dock
            self.dock_pending = np.zeros(g.size)
        if shorten:
            s_rate = self._shortening_rate(src)                     # m/yr per overriding trench node
            moved = margin_shortening(near, d, ok & cont, s_rate[np.clip(near, 0, None)], w['thick']+out, g.area, dt,
                                      p['subduction_d0_cont_km']*1e3, p['subduction_sigma_cont_km']*1e3,
                                      p['foreland_d0_km']*1e3, p['foreland_sigma_km']*1e3,
                                      match=p['shortening_match'])
            out = out+moved
            if p['belt_lithosphere'] == 'shorten':
                self.belt_gain = np.where(self.belt, np.maximum(moved, 0.), 0.)
            self.erosion_totals['margin_shortening_moved_m3'] = self.erosion_totals.get('margin_shortening_moved_m3', 0.)+float(
                (np.maximum(moved, 0)*g.area).sum())
            self._record_upper_plate(src)
            if p['backarc_extension']:
                out = out+self._backarc_extension(near, d, w['thick']+out)
        elif p['backarc_extension']:
            raise ValueError('backarc_extension needs margin_shortening (the same upper-plate rate)')
        return out


    def _record_upper_plate(self, src):
        """Output only: overriding trench length (m, summed over steps) in all, on continental crust, and of that in
        shortening / extension by the rate u - rollback (and beyond +-5 mm/yr), for Earth's share of transects in each
        regime (Lallemand et al. 2008)."""
        g, w = self.g, self.world
        rate, front, ell = self._upper_plate_rate(src)
        ET = self.erosion_totals
        has = front >= 0
        c = has & w['cont']
        for key, m in (('overriding_trench_m', has), ('overriding_cont_trench_m', c),
                       ('overriding_cont_shortening_m', c & (rate > 0)), ('overriding_cont_extending_m', c & (rate < 0)),
                       ('overriding_cont_shortening5_m', c & (rate > 5e-3)), ('overriding_cont_extending5_m', c & (rate < -5e-3)),
                       ('overriding_ocean_extending_m', has & ~w['cont'] & (rate < 0))):
            ET[key] = ET.get(key, 0.)+float(ell[m].sum())


    def _backarc_extension(self, near, d, thick):
        """backarc_extension (see PARAMS): thinning (m, <= 0) of the back-arc bands behind the trench segments whose
        arc-forearc block steps trenchward this step; the steps are kept in self.backarc_steps for _backarc_strips."""
        g, p, w = self.g, self.p, self.world
        out = np.zeros(g.size)
        self.backarc_steps = []
        src = (self.btype == SUBDUCTION) & ~self.trench
        rate, front, ell = self._upper_plate_rate(src)
        e = np.where((front >= 0) & w['cont'] & (rate < 0), -rate, 0.)
        I = np.flatnonzero(e > 0)
        if I.size == 0:
            return out
        if not hasattr(self, 'backarc_rng'):
            self.backarc_rng = np.random.default_rng([self.seed, 10])
        a = g.area
        dt = p['dt_myr']*1e6
        F = front[I]
        P = e[I]*dt*ell[I]/a[F]
        ET = self.erosion_totals
        ET['backarc_expected_m2'] = ET.get('backarc_expected_m2', 0.)+float((e[I]*dt*ell[I]).sum())
        ET['backarc_capped_steps'] = ET.get('backarc_capped_steps', 0)+int((P > 1).sum())
        go = self.backarc_rng.random(I.size) < P
        I, F = I[go], F[go]
        if I.size == 0:
            return out
        _, first = np.unique(F, return_index=True)           # one step per front node
        keep = np.zeros(I.size, bool); keep[first] = True
        I, F = I[keep], F[keep]
        idx = self._claim_idx
        slot = self.owner_slot
        def free(i, f):                                  # the front's raster node is new to the plate and touches it
            m = self.plates[slot[i]].mask
            j = idx[slot[i], f]
            nj = g.nbr[j]
            return (not m[j]) and bool(m[nj[nj >= 0]].any())
        ok = np.array([free(i, f) for i, f in zip(I, F)], bool)
        I, F = I[ok], F[ok]
        if I.size == 0:
            return out
        # the back-arc band of each stepping segment: continental nodes of the same plate whose nearest overriding
        # trench node is the segment's, weighted by Earth's profile at their distance from the trench
        km = np.asarray(BACKARC_PROFILE_KM)*1e3
        b = np.searchsorted(km, d, side='right')-1
        inb = (b >= 0) & (b < len(BACKARC_PROFILE_W))
        k = np.where(inb, np.asarray(BACKARC_PROFILE_W)[np.clip(b, 0, len(BACKARC_PROFILE_W)-1)], 0.)
        nr = np.clip(near, 0, None)
        Tb = p['breakup_thick_m']
        zone = (near >= 0) & inb & w['cont'] & (self.owner == self.owner[nr]) & (thick > Tb)
        sel = np.zeros(g.size, bool); sel[I] = True
        zone &= sel[nr]
        K = np.bincount(near[zone], weights=(a*k)[zone], minlength=g.size)
        has = K[I] > 0
        I, F = I[has], F[has]
        if I.size == 0:
            return out
        af = np.zeros(g.size); af[I] = a[F]
        zone &= af[nr] > 0
        c = np.where(zone, af[nr]*k/np.maximum(K[nr], 1e-300), 0.)
        out = np.where(zone, np.maximum(-thick*(1.-np.exp(-c)), Tb-thick), 0.)
        V = np.bincount(near[zone], weights=(-out*a)[zone], minlength=g.size)
        if p['lithosphere'] == 'thermal_age':        # C: the stretched lithosphere thins with the crust (as rifts)
            th = zone & (out < 0)
            ratio = np.clip((thick[th]+out[th])/thick[th], 0., 1.)
            w['therm'][th] = self.t-stretched_age(self.t-w['therm'][th], ratio, p)
            self.therm_changed |= th
        self.backarc_steps = [(int(i), int(f), float(V[i])) for i, f in zip(I, F)]
        ET['backarc_steps'] = ET.get('backarc_steps', 0)+int(I.size)
        ET['backarc_stretched_m3'] = ET.get('backarc_stretched_m3', 0.)+float(V[I].sum())
        ET['backarc_area_m2'] = ET.get('backarc_area_m2', 0.)+float(a[F].sum())
        return out


    def _backarc_strips(self):
        """backarc_extension: after the rasters' update, each stepping segment's plate claims its front node with the
        band's stretched volume as continental crust (the translated arc-forearc block; it keeps the edge node's crust,
        orogeny and thermal ages; no sediment). Ledgers 'backarc' (volume) and area 'backarc'."""
        steps = getattr(self, 'backarc_steps', [])
        self.backarc_steps = []
        if not steps:
            return
        g, w = self.g, self.world
        a = g.area
        idx = self._claim_idx
        for i, f, V in steps:
            k = self.owner_slot[i]
            pl = self.plates[k]
            j = idx[k, f]
            T = V/a[j]
            if pl.mask[j] and pl.cont[j]:            # covered since the step was drawn: the volume joins that node
                pl.thick[j] += np.float32(T)
            else:
                pl.mask[j] = True
                pl.cont[j] = True
                pl.thick[j] = T
                pl.born[j] = w['born'][i]
                pl.orog[j] = w['orog'][i]
                pl.therm[j] = w['therm'][i]
                pl.sed[j] = 0.
                pl.under[j] = 0
                pl.slab[j] = 0
                self.area_ledger['backarc'] = self.area_ledger.get('backarc', 0.)+float(a[j])
            self.ledger['backarc'] = self.ledger.get('backarc', 0.)+float(np.float32(T))*float(a[j])


    def _thermal_resets(self, collided):
        """C: continental world nodes whose lithosphere is reset to thermal age 0 this step: collision contacts (crust
        consumed under continental crust, or a collision boundary closing faster than boundary_min_mm_yr; not the belt
        clock, which resets on thickness alone) and the continental arc belt over a converging trench (with
        belt_lithosphere 'shorten', the belt nodes that margin shortening did not thicken this step)."""
        g, p, w = self.g, self.p, self.world
        cont = w['cont']
        belt = cont & getattr(self, 'belt', np.zeros(g.size, bool))
        if p.get('belt_lithosphere', 'reset') == 'shorten':  # shortened belt nodes thickened instead (sim.step)
            belt = belt & ~(getattr(self, 'belt_gain', np.zeros(g.size)) > 0)
        if p.get('collision_lithosphere', 'reset') == 'thicken':    # collisions thicken the lithosphere instead
            return belt
        out = cont & collided
        if hasattr(self, 'btype') and hasattr(self, 'conv_rate'):
            out |= cont & (self.btype == COLLISION) & (self.conv_rate*1e3 > p['boundary_min_mm_yr'])
        return out | belt


    def _orogen_marks(self):
        """World nodes (continental crust) besides collision contacts whose orogeny clock is reset this step
        (orogen_rule; collision contacts are always marked)."""
        g, p, w = self.g, self.p, self.world
        rule = p['orogen_rule']
        none = np.zeros(g.size, bool)
        if rule == 'collision':
            return none
        cont = w['cont']
        if rule == 'belt':
            return cont & (getattr(self, 'belt', none) | (w['thick'] > p['collapse_thick_m']))
        if rule == 'strain':
            if not hasattr(self, 'tectonic'):
                return none
            rate = np.abs(self.tectonic)/np.maximum(w['thick'], 1e3)/(p['dt_myr']*1e6)
            return cont & (rate > p['orogen_strain_rate'])
        raise ValueError(f"orogen_rule {rule!r}: 'collision', 'belt' or 'strain'")


    def _shortening_rate(self, src):
        """Upper-plate shortening rate (m/yr) at every overriding trench node: max(0, u - rollback), u the overriding
        plate's speed toward the trench in the mantle (world) frame (Oncken et al. 2006). The trench direction is the
        mean direction to the node's neighbours on other plates."""
        g, p = self.g, self.p
        out = np.zeros(g.size)
        idx = np.flatnonzero(src)
        nb = g.nbr[idx]
        other = (nb >= 0) & (self.owner[np.maximum(nb, 0)] != self.owner[idx][:, None])
        u, n = self._trenchward_speed(idx, nb, other)
        out[idx] = np.clip(u-self._rollback(idx, nb, other, n), 0, None)*other.any(1)
        return out


    def _rollback(self, idx, nb, other, n):
        """The trench's rollback (m/yr, oceanward) at overriding trench nodes idx (rollback_rule, see PARAMS)."""
        p = self.p
        rule = p['rollback_rule']
        if rule == 'fixed':
            return p['rollback_mm_yr']*1e-3
        if rule != 'lallemand':
            raise ValueError(f"rollback_rule {rule!r}: 'fixed' or 'lallemand'")
        g, w = self.g, self.world
        nbs = np.maximum(nb, 0)
        oc = other & ~w['cont'][nbs]
        use = np.where(oc.any(1)[:, None], oc, other)       # the subducting plate: oceanic neighbours, else any
        om = {pl.id: pl.omega(g.radius) for pl in self.plates}
        x = g.xyz[idx]
        vs = np.zeros(idx.size)
        for c in range(nb.shape[1]):
            wj = np.stack([om.get(int(k), np.zeros(3)) for k in self.owner[nbs[:, c]]]) if idx.size else np.zeros((0, 3))
            vs += np.where(use[:, c], -np.einsum('ij,ij->i', np.cross(wj, x)*g.radius, n), 0.)
        vs /= np.maximum(use.sum(1), 1)                       # m/yr toward the upper plate
        return (p['rollback_sb04_mm_yr']*1e-3-p['rollback_sb04_slope']*vs)


    def _trenchward_speed(self, idx, nb, other):
        """u (m/yr): the overriding plate's speed toward the trench at nodes idx (world frame), and the unit trench
        normal n (mean direction to the neighbours on other plates)."""
        g = self.g
        x = g.xyz[idx]
        disp = np.where(other[..., None], g.xyz[np.maximum(nb, 0)]-x[:, None], 0.).sum(1)
        disp -= np.einsum('ij,ij->i', disp, x)[:, None]*x
        n = disp/np.maximum(np.linalg.norm(disp, axis=1, keepdims=True), 1e-15)
        om = {pl.id: pl.omega(g.radius) for pl in self.plates}
        w = np.stack([om.get(int(k), np.zeros(3)) for k in self.owner[idx]]) if idx.size else np.zeros((0, 3))
        u = np.einsum('ij,ij->i', np.cross(w, x)*g.radius, n)             # m/yr toward the trench
        return u, n


    def _upper_plate_rate(self, src):
        """The signed upper-plate rate u - rollback (m/yr; > 0 shortening, < 0 extension; Heuret & Lallemand 2005) at
        every overriding trench node with an oceanic neighbour on another plate, that neighbour (the node most nearly
        along the trench normal: the front the arc-forearc block steps onto; -1 elsewhere), and the trench length per
        node (BOUNDARY_NODE_H h)."""
        from planet_tectonics.motion import BOUNDARY_NODE_H
        g, p, w = self.g, self.p, self.world
        rate = np.zeros(g.size)
        front = np.full(g.size, -1, np.int64)
        idx = np.flatnonzero(src)
        nb = g.nbr[idx]
        nbs = np.maximum(nb, 0)
        other = (nb >= 0) & (self.owner[nbs] != self.owner[idx][:, None])
        u, n = self._trenchward_speed(idx, nb, other)
        rate[idx] = (u-self._rollback(idx, nb, other, n))*other.any(1)
        ocean = other & ~w['cont'][nbs]
        along = np.einsum('ikj,ij->ik', g.xyz[nbs]-g.xyz[idx][:, None], n)
        best = np.argmax(np.where(ocean, along, -np.inf), 1)
        has = ocean.any(1)
        front[idx[has]] = nbs[has, best[has]]
        ell = BOUNDARY_NODE_H*np.sqrt(g.area)
        return rate, front, ell
