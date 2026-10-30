"""PlateSim, start: the supercontinent and its plates, start poles."""
import numpy as np

from planet_grid import fields
from planet_landscape import graph_landscape as gl

from planet_tectonics.crust_kernels import quantile_map
from planet_tectonics.graphops import components, nearest_label
from planet_tectonics.params import EARTH_STABLE_CRUST_KM, RIDGE
from planet_tectonics.plate import Plate


class StartMixin:
    """PlateSim, start: the supercontinent and its plates, start poles."""

    # ------------------------------------------------------------------ start
    def _new_plate(self, axis, speed):
        pl = Plate(self.next_id, self.g.size, axis, speed)
        self.next_id += 1
        self.plates.append(pl)
        return pl


    def _init_supercontinent(self, n_plates, n_cont_plates, frac):
        g, rng, p = self.g, self.rng, self.p
        xyz = g.xyz
        s = rng.normal(size=3); s /= np.linalg.norm(s)
        self.sc_centre = s
        ang = np.arccos(np.clip(xyz@s, -1, 1))
        score = ang+.35*fields.random_field(xyz, 12, 2.5, rng)
        cont = score <= np.quantile(score, frac)
        lab = components(cont, g.nbr)
        sizes = np.bincount(lab[lab >= 0], weights=g.area[lab >= 0])
        cont = lab == np.argmax(sizes)
        ocean_lab = components(~cont, g.nbr)            # fill small enclosed seas
        osz = np.bincount(ocean_lab[ocean_lab >= 0], weights=g.area[ocean_lab >= 0])
        small = np.flatnonzero(osz < .02*g.area.sum())
        cont |= np.isin(ocean_lab, small)

        # cratons (Voronoi cells inside the continent) welded by old orogens
        cids = np.flatnonzero(cont)
        n_cr = int(rng.integers(7, 12))
        seeds = np.full(g.size, -1, np.int64)
        seeds[rng.choice(cids, n_cr, replace=False)] = np.arange(n_cr)
        clab, _ = nearest_label(np.where(cont, seeds, -2), g.nbr, np.where(cont[g.nbr] & cont[:, None], g.dist, 1e30))
        boundary = self._label_boundary(cont, clab)
        d_suture = gl.distance_to(boundary, g.nbr, g.dist)
        d_coast = gl.distance_to(~cont, g.nbr, g.dist)
        craton_age = rng.uniform(1500, 3500, n_cr)
        pair_age = rng.uniform(300, 1500, (n_cr, n_cr)); pair_age = np.triu(pair_age)+np.triu(pair_age, 1).T
        thick = (p['cont_ref_thick_m']-500.+1500.*fields.random_field(xyz, 24, 2., rng)
                 + 4000.*np.exp(-d_suture/150e3)+2000.*np.clip((d_suture-300e3)/300e3, 0, 1))
        if p['start_crust'] == 'earth':                 # Earth's stable crust distribution, the field's order kept
            thick[cont] = quantile_map(thick[cont], g.area[cont], target=tuple(1e3*k for k in EARTH_STABLE_CRUST_KM))
        elif p['start_crust'] != 'assumed':
            raise ValueError(f"start_crust {p['start_crust']!r}: 'assumed' or 'earth'")
        thick = thick-(p['cont_ref_thick_m']-28000.)*np.clip(1-d_coast/200e3, 0, 1)   # thinned rim at the margins
        born_c = self.t-craton_age[np.clip(clab, 0, None)]
        # old orogen age: from the two cratons it welds (nearest other-craton label across the suture)
        other = self._other_label(cont, clab)
        orog_c = np.where(d_suture < 250e3, self.t-pair_age[np.clip(clab, 0, None), np.clip(other, 0, None)], np.nan)

        # plates: Voronoi with wiggled distances; n_cont_plates seeds inside the continent
        oids = np.flatnonzero(~cont)
        pseeds = np.concatenate([rng.choice(cids, n_cont_plates, replace=False), rng.choice(oids, n_plates-n_cont_plates, replace=False)])
        lab0 = np.full(g.size, -1, np.int64); lab0[pseeds] = np.arange(n_plates)
        wig = 1+.35*fields.random_field(xyz, 10, 2., rng)
        w = np.clip(.5*(wig[:, None]+wig[np.clip(g.nbr, 0, None)]), .3, None)
        if p['start_plates'] == 'voronoi':
            plab, _ = nearest_label(lab0, g.nbr, g.dist*w)
        elif p['start_plates'] == 'coast':
            nbc = np.clip(g.nbr, 0, None)
            same = (g.nbr >= 0) & (cont[nbc] == cont[:, None])
            dw = np.where(same, g.dist*w, np.inf)
            plab, _ = nearest_label(lab0, g.nbr, dw)
            while (plab < 0).any():                      # a domain component without a seed: its own plate
                un = plab < 0
                comp = components(un, g.nbr)
                first = comp == comp[np.flatnonzero(un)[0]]
                dc = gl.distance_to(~first, g.nbr, g.dist)
                lab0[int(np.argmax(np.where(first, dc, -1)))] = n_plates
                n_plates += 1
                plab, _ = nearest_label(lab0, g.nbr, dw)
        else:
            raise ValueError(f"start_plates {p['start_plates']!r}: 'voronoi' or 'coast'")

        # ocean ages: distance from initial ridges (ocean-ocean plate boundaries), half rate 3 cm/yr, cap 180 Myr
        ridge = (~cont) & self._label_boundary(~cont, plab)
        d_ridge = gl.distance_to(ridge, g.nbr, g.dist)
        age_o = np.clip(d_ridge/30e3*(1+.2*fields.random_field(xyz, 16, 2., rng)), 0, 180)

        for k in range(n_plates):
            own = plab == k
            c = (g.area[own, None]*xyz[own]).sum(0); c /= np.linalg.norm(c)
            if k < n_cont_plates:
                axis = np.cross(s, c)+.3*rng.normal(size=3)*np.linalg.norm(np.cross(s, c))
                speed = rng.uniform(2., 4.)
            else:
                axis = rng.normal(size=3)
                speed = rng.uniform(3., 7.)
            pl = self._new_plate(axis, speed)
            pl.mask[:] = own
            pl.cont[:] = own & cont
            pl.thick[:] = np.where(cont, thick, p['ocean_thick_m'])
            pl.born[:] = np.where(cont, born_c, self.t-age_o)
            pl.orog[:] = np.where(cont, orog_c, np.nan)


    def _start_poles_from_forces(self):
        """start_poles 'forces': oceanic start plates move as the force law sets them on the start's geometry, with
        coast nodes on the ocean side as trench and ocean-ocean boundary nodes as ridge."""
        g = self.g
        nb = g.nbr; ok = nb >= 0; nbc = np.where(ok, nb, 0)
        c = self.world['cont']
        other = ok & (self.owner[nbc] != self.owner[:, None])
        trench = (~c) & np.any(other & c[nbc], 1)
        ridge = (~c) & np.any(other & ~c[nbc], 1) & ~trench
        self.trench, self.btype = trench, np.where(ridge, RIDGE, 0).astype(np.int8)
        for pl in self.plates:
            own = self.owner == pl.id
            if own.any() and not (own & c).any():
                pl.speed, pl.axis = self._driving(own, pl.speed, pl.axis, relax=1.)
        del self.trench, self.btype                      # typed by motion from the first step on
