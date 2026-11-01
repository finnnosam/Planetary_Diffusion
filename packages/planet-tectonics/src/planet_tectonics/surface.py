"""PlateSim, surface: coarse erosion on the composed surface and the isostatic and flexural elevation."""
import numpy as np

from planet_landscape import graph_landscape as gl
from planet_tectonics import shelf as shelf_mod

from planet_tectonics.isostasy import gdh1_depth, sea_level_for_volume, sediment_column


class SurfaceMixin:
    """PlateSim, surface: coarse erosion on the composed surface and the isostatic and flexural elevation."""

    def _erode(self, thick):
        """Coarse erosion over one plate step: the landscape model's surface processes (graph_landscape.surface_step)
        on the flexural-isostatic elevation, in erosion_substeps steps. Returns the change of crust thickness
        (continental: erosion < 0, and without continental_sediment sediment deposited on continental crust > 0;
        oceanic, with ocean_crust_erosion: erosion below the sediment, <= 0) and of sediment (ocean floor, and with
        continental_sediment continental crust too: deposits > 0, erosion of the sediment first). Their sum times area
        is 0 (surface_step is closed)."""
        g, p, w = self.g, self.p, self.world
        z0 = self.elevation(thick, w['sed'])
        self.z_last = z0                                 # dyn_rule 'slabs' water-loads the flow where this is < 0
        lp = {**gl.DEFAULTS, 'K': gl.DEFAULTS['K']*p['erosion_K_factor'], 'Sc': p['erosion_Sc'],
              'Kd': gl.DEFAULTS['Kd']*p['erosion_Kd_factor'], 'open_edges': False, 'theta': p['erosion_theta']}
        since = np.where(np.isnan(w['orog']), np.inf, self.t-w['orog'])
        erod = p['stable_erodibility']+(1-p['stable_erodibility'])*np.exp(-np.clip(since, 0, None)/p['orogen_erodibility_tau_myr'])
        K = lp['K']*np.where(w['cont'], erod, 1.)
        seed = z0 < -1500.
        dt = p['dt_myr']*1e6/p['erosion_substeps']
        z = z0.copy()
        fallback = np.ones(g.size, bool)
        shelf_on = p['shelf_kappa_m2_yr'] > 0
        keep = p['keep_deep_sediment']
        deep = np.zeros(g.size)
        for _ in range(p['erosion_substeps']):
            r = gl.surface_step(g, z, seed, None, fallback, K, lp, dt, self.h_min, self.bound, keep_deep=keep)
            z = r['z']
            if keep and r['deep_m3'].any():          # turbidites: fill the sea-floor lows from the bottom up
                deep += r['deep_m3']
                z = gl.fill_and_spill(z, r['deep_m3'], g.nbr, g.area)
            if shelf_on:                             # wave planation and coastal retreat (planet_tectonics.shelf)
                ocean = gl.connected_ocean(z, seed, g.nbr)
                avail = np.where(w['cont'], 1e9, np.clip(w['sed']+z-z0, 0, None))
                z, _, off = shelf_mod.shelf_step(g, z, ocean, avail, p['shelf_kappa_m2_yr'], p['shelf_wave_base_m'],
                                                 p['cliff_retreat_m_yr'], dt)
                self.erosion_totals['shelf_off_land_m3'] = self.erosion_totals.get('shelf_off_land_m3', 0.)+off
        self.erosion_totals['deep_sea_sediment_m3'] = self.erosion_totals.get('deep_sea_sediment_m3', 0.)+float(deep.sum())
        dz = z-z0
        cont = w['cont']
        a = g.area
        # surface change on continental nodes (sediment and crust together): the totals keep their meaning
        dz_cont = np.where(cont, dz, 0.)
        self.erosion_totals['continental_eroded_m3'] += float((np.minimum(dz_cont, 0)*a).sum())
        self.erosion_totals['continental_deposited_m3'] += float((np.maximum(dz_cont, 0)*a).sum())
        # the part deposited on continental crust under the sea (shelves, drowned margins): delivered to the sea in
        # Earth's sense (river delivery, 4-8 km^3/yr), so scripts.sweep_plates counts it as export
        self.erosion_totals['continental_deposited_submerged_m3'] = self.erosion_totals.get(
            'continental_deposited_submerged_m3', 0.)+float((np.maximum(dz_cont, 0)*a)[z0 <= 0].sum())
        if p['continental_sediment']:                # deposits are sediment everywhere; erosion takes sediment first
            dsed = np.maximum(dz, -w['sed'])
            e_cont = np.where(cont, dz-dsed, 0.)
            self.erosion_totals['continental_sediment_m3'] = self.erosion_totals.get('continental_sediment_m3', 0.)+float(
                (dsed*a)[cont].sum())
            dsed_oc = np.where(cont, 0., dsed)
        else:                                        # deposits on continental crust become crust
            e_cont = dz_cont
            dsed = dsed_oc = np.where(cont, 0., np.maximum(dz, -w['sed']))
        self.erosion_totals['continental_crust_eroded_m3'] = self.erosion_totals.get('continental_crust_eroded_m3', 0.)+float(
            (np.minimum(e_cont, 0)*a).sum())
        self.erosion_totals['ocean_sediment_m3'] += float((dsed_oc*a).sum())
        if p['ocean_crust_erosion']:                 # cut below the sediment: from the oceanic crust (closes the ledger)
            e_oc = np.where(cont, 0., np.maximum(np.minimum(dz+w['sed'], 0.), -thick))
            self.erosion_totals['ocean_crust_eroded_m3'] = self.erosion_totals.get('ocean_crust_eroded_m3', 0.)+float(
                (e_oc*a).sum())
            self.erosion_totals['ocean_crust_eroded_subaerial_m3'] = self.erosion_totals.get(
                'ocean_crust_eroded_subaerial_m3', 0.)+float((e_oc*a)[z0 > 0].sum())
            normal = (z0 > 0) & (thick < p['ocean_thick_m']+1000.)    # the old leak: normal oceanic crust above sea
            self.erosion_totals['ocean_crust_eroded_subaerial_normal_m3'] = self.erosion_totals.get(
                'ocean_crust_eroded_subaerial_normal_m3', 0.)+float((e_oc*a)[normal].sum())
            e_cont = e_cont+e_oc
        return e_cont, dsed


    def elevation(self, thick=None, sed=None):
        """Elevation (m) with flexural isostasy on an elastic plate (Te).

        Continental crust: its thickness anomaly (relative to cont_ref_thick_m, which stands at cont_ref_elev_m) is a
        load of density rho_crust against the mantle. Ocean floor: thermal subsidence with age (GDH1); its crust
        anomaly (arcs, plateaus) and sediment are submarine loads of density (rho - rho_water) against
        (rho_mantle - rho_water), so a basin filled to sea level holds ~2.5x its depth of sediment, as on Earth.
        Long wavelengths reduce to Airy isostasy; valleys and single belts are partly held up by the plate.
        A thinned continental column below sea level is water-loaded (its own Airy depth only, not its neighbours').
        With thinning_local_isostasy, continental crust thinner than cont_ref_thick_m stands at that local Airy depth
        (plus the deflection from other loads) and is no negative load on the plate. With continental_sediment, the
        sediment on continental crust adds its thickness and loads the plate (rho_sediment; against water where its
        top is below sea level). With sediment_compaction the sediment's height and mass are those of a compacting
        column (sediment_column): bulk height on the surface, solid mass and pore water in the load.
        """
        p, w = self.p, self.world
        thick = w['thick'] if thick is None else thick
        sed = w['sed'] if sed is None else sed
        cont = w['cont']
        rc, rm, rw = p['rho_crust'], p['rho_mantle'], p['rho_water']
        sh, sm_dry, sm_wet = sediment_column(sed, cont, p)       # height, dry mass, submerged excess mass (per area)
        anomaly = np.where(cont, thick-p['cont_ref_thick_m'], thick-p['ocean_thick_m'])
        ocean_age = np.where(cont, 0., self.t-w['born'])
        s0 = getattr(self, 'sea_level', 0.)                     # previous sea level on the datum (loads use it)
        above = (-gdh1_depth(ocean_age)+anomaly+sh) > s0        # oceanic crust built above sea level (arcs, fills)
        # thinning_local_isostasy: the thinned part (negative anomaly) of continental crust stands at its local Airy
        # level and is no load on the plate; the thickened part loads it as before
        thin = np.where(cont, np.minimum(anomaly, 0.), 0.) if p['thinning_local_isostasy'] else 0.
        load = anomaly-thin
        own_airy = p['cont_ref_elev_m']+anomaly*(rm-rc)/rm-s0
        # continental_sediment: sediment on continental crust is a load like ocean-floor sediment, dry where its top
        # (on the crust's own water-loaded Airy level, the ocean rule's counterpart) stands above sea level
        zero = np.zeros_like(anomaly)
        cs = np.where(cont, sh, 0.) if p['continental_sediment'] else zero
        cs_dry = np.where(cont, sm_dry, 0.) if p['continental_sediment'] else zero
        cs_wet = np.where(cont, sm_wet, 0.) if p['continental_sediment'] else zero
        c_dry = cont & (np.where(own_airy < 0, own_airy*rm/(rm-rw), own_airy)+cs > 0)
        dry = np.where(cont, rc*load+np.where(c_dry, cs_dry, 0.), np.where(above, rc*anomaly+sm_dry, 0.))
        wet = np.where(cont, np.where(c_dry, 0., cs_wet), np.where(above, 0., (rc-rw)*anomaly+sm_wet))
        defl = self.flex_dry.deflection(dry)+self.flex_wet.deflection(wet)
        z = np.where(cont, p['cont_ref_elev_m']+load+thin*(rm-rc)/rm+cs, -gdh1_depth(ocean_age)+anomaly+sh)-defl
        z = z+np.where(cont & (own_airy < 0), own_airy*rw/(rm-rw), 0.)
        z = z+self.dyn                                    # dynamic topography: non-load surface motion
        if p['water_volume_km3'] is None:
            return z
        self.sea_level = sea_level_for_volume(z, self.g.area, p['water_volume_km3']*1e9)
        return z-self.sea_level
