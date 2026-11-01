"""Mantle-frame forcing for the plate layer (docs/planet-layer-plan.md, sect. 3): dynamic topography and hotspots.

The world frame of PlateSim is the mantle frame, so plates drift over these fields and each crust parcel feels them
as it passes (which is what the landscape stage's parcel histories need).

Dynamic topography: a long-wavelength field (spherical-harmonic degrees 1-12, red spectrum) that evolves by blending
independent random fields every `dt_knot_myr`; standard deviation `dyn_std_m`. It is *non-load* vertical motion: it
moves the surface without loading the plate (residual topography on Earth is ~ +-1 km, changing at ~10-100 m/Myr;
Hoggard et al. 2016). Hotspots add a swell (part of the dynamic topography) and volcanic crust to whatever plate
passes over them; a new plume starts with a large igneous province (a burst of crust over a wide disc).

With dyn_rule 'slabs' (D, journal 66) the knot field is replaced by the free-slip flow of the planet's own slab record
(planet_tectonics.mantle_flow); the knots are still drawn, so the hotspots take the same random numbers in both rules.
"""
import numpy as np

from planet_grid import fields

PARAMS = {
    'dyn_std_m': 400., 'dyn_lmax': 12, 'dyn_slope': 2.5, 'dt_knot_myr': 40.,
    'hotspots_initial': 10, 'hotspot_birth_per_myr': .05,
    'hotspot_life_myr': (60., 180.),
    'hotspot_flux_km3_myr': (2000., 15000.),      # steady volcanic output (Hawaii ~1e6 km3 over ~70 Myr)
    'hotspot_radius_km': 60., 'hotspot_swell_m': 800., 'hotspot_swell_radius_km': 500.,
    'lip_volume_km3': (5e5, 4e6), 'lip_radius_km': 600., 'lip_duration_myr': 2.,
    'dyn_rule': 'random', 'rho_mantle': 3300., 'rho_water': 1030., 'threads': None,
    'slab_dip_deg': 58., 'lm_sink_mm_yr': 12., 'flow_z_cut_km': 300., 'flow_lmax': 31, 'viscosity_file': None,
    'slab_alpha': 'thermo',
}


class Mantle:
    def __init__(self, graph, rng, t0, t_end=0., params=None):
        self.g, self.rng = graph, rng
        self.p = {**PARAMS, **(params or {})}
        p = self.p
        self.knots_t = np.arange(t0-p['dt_knot_myr'], t_end+2*p['dt_knot_myr'], p['dt_knot_myr'])
        self.knots = [p['dyn_std_m']*fields.random_field(graph.xyz, p['dyn_lmax'], p['dyn_slope'], rng)
                      for _ in self.knots_t]
        if p['dyn_rule'] == 'slabs':
            from planet_tectonics.mantle_flow import MantleFlow
            self.flow = MantleFlow(graph, p)
        elif p['dyn_rule'] == 'random':
            self.flow = None
        else:
            raise ValueError(f"dyn_rule {p['dyn_rule']!r}: 'random' or 'slabs'")
        self.hotspots = []
        for _ in range(p['hotspots_initial']):
            self._new_hotspot(t0, lip=False, age=rng.uniform(0, p['hotspot_life_myr'][0]))

    def _new_hotspot(self, t, lip=True, age=0.):
        p, rng = self.p, self.rng
        v = rng.normal(size=3); v /= np.linalg.norm(v)
        life = rng.uniform(*p['hotspot_life_myr'])
        self.hotspots.append({'pos': v, 'start': t-age, 'end': t-age+life,
                              'flux_m3_yr': np.exp(rng.uniform(*np.log(p['hotspot_flux_km3_myr'])))*1e9/1e6,
                              'lip_m3': np.exp(rng.uniform(*np.log(p['lip_volume_km3'])))*1e9 if lip else 0.})

    def advance(self, t, dt):
        """Births and deaths of plumes over (t - dt, t]."""
        p = self.p
        for _ in range(self.rng.poisson(p['hotspot_birth_per_myr']*dt)):
            self._new_hotspot(t)
        self.hotspots = [h for h in self.hotspots if h['end'] > t]

    def _disc(self, pos, radius_m):
        ang = np.arccos(np.clip(self.g.xyz@pos, -1, 1))*self.g.radius
        return np.exp(-.5*(ang/radius_m)**2)

    def subduct(self, t, dt, xyz, dirn, mass_kg, speed_m_yr):
        """dyn_rule 'slabs': move the record over (t - dt, t], add this step's parcels (at the trench, z = 0), merge
        the lower mantle and solve the flow."""
        rec = self.flow.record
        rec.advance(dt, self.g.radius)
        rec.deposit(xyz, dirn, mass_kg, speed_m_yr, t)
        rec.merge(self.g, self.flow.z_cut)
        self.flow.update()

    def swells(self, t):
        out = np.zeros(self.g.size)
        for h in self.hotspots:
            if h['start'] <= t < h['end']:
                out = out+self.p['hotspot_swell_m']*self._disc(h['pos'], self.p['hotspot_swell_radius_km']*1e3)
        return out

    def topography(self, t, wet=None):
        """Dynamic topography (m) at time t: blended random fields (dyn_rule 'random') or the slab flow's air-loaded
        field, x rho_m / (rho_m - rho_w) on the nodes in wet (dyn_rule 'slabs'), plus hotspot swells."""
        if self.flow is not None:
            p = self.p
            f = 1. if wet is None else np.where(wet, p['rho_mantle']/(p['rho_mantle']-p['rho_water']), 1.)
            return self.flow.h_air*f+self.swells(t)
        k = np.clip(np.searchsorted(self.knots_t, t)-1, 0, len(self.knots_t)-2)
        f = np.clip((t-self.knots_t[k])/(self.knots_t[k+1]-self.knots_t[k]), 0, 1)
        f = f*f*(3-2*f)
        out = (1-f)*self.knots[k]+f*self.knots[k+1]
        for h in self.hotspots:
            if h['start'] <= t < h['end']:
                out = out+self.p['hotspot_swell_m']*self._disc(h['pos'], self.p['hotspot_swell_radius_km']*1e3)
        return out

    def volcanic(self, t, dt):
        """Volcanic crust added over (t - dt, t] (m of crust per node): steady plume output and LIP bursts."""
        p, a = self.p, self.g.area
        out = np.zeros(self.g.size)
        for h in self.hotspots:
            if not (h['start'] < t and h['end'] > t-dt):
                continue
            w = self._disc(h['pos'], p['hotspot_radius_km']*1e3)
            out += h['flux_m3_yr']*dt*1e6*w/(w*a).sum()
            if h['lip_m3'] > 0 and t-dt < h['start']+p['lip_duration_myr']:
                wl = self._disc(h['pos'], p['lip_radius_km']*1e3)
                out += h['lip_m3']*min(dt, p['lip_duration_myr'])/p['lip_duration_myr']*wl/(wl*a).sum()
        return out
