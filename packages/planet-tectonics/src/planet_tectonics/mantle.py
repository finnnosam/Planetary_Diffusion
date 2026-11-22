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
With flow_power 'measured' (ninety-ninth pass, journal 132) the flow keeps its pattern and takes Earth's measured power
per degree (VD2020_POWER_KM2 below, air-loaded): the measured field placed by the plates' own slab history. The record
starts empty, so the start's mantle is dyn_rule 'measured''s own field (the same stream, so the start is the control's)
and it hands over to the placed field at constant variance over one knot interval (dt_knot_myr), the measured field's
own decorrelation step: no new number. (Without it the start had no field, and start_cover 'accommodation''s regression
on a bare level nearly collinear with the crust put Earth's thickest cover on thick crust: a 9.5 km start node.)

With flow_placed 'earth_corr' (pass 101, journal 136) the placed flow keeps only the share of each degree's measured
power that the same flow on Earth's own subduction history accounts for on Earth: max(r_l, 0)^2, r_l the degree
correlation of Earth's measured residual topography (VD2020's posterior mean) with E3 (Muller et al. 2019's slab
history through D's record and kernels, journal 70; EARTH_SLAB_CORR below). The rest of the power, (1 - share), is a
plate-blind field like dyn_rule 'measured''s, drawn from the same stream after the start's field and blended between
knots at constant variance. Earth's slabs account for l 1-2 (r 1.00, 0.88), little of l 3 (0.32) and nothing above.

With dyn_rule 'measured' (fifty-eighth pass, journal 80) the knots are random fields with Earth's measured residual-
topography spectrum per degree (VD2020_POWER_KM2: Valentine & Davies 2020), drawn from their own stream (field_rng)
so the hotspots again take the same random numbers. Nothing is fitted: no amplitude and no slope. The measurement is
water-loaded; the field is held air-loaded (x 2,300 / 3,300, the measurement's own load densities) and water-loaded
on the nodes under the sea as D's flow is (x rho_m / (rho_m - rho_w)). Between knots the fields blend at constant
variance (cos / sin of an angle that runs smoothly from 0 to pi/2), so the amplitude does not pulse.
"""
import numpy as np

from planet_grid import fields

# Earth's residual topography, expected power per degree l = 1 .. 30 (km^2; orthonormal real harmonics, sum over m;
# water-loaded with Delta rho 2,300 kg/m^3): Valentine & Davies (2020), G-cubed 21, e2020GC009240, the posterior of their
# Gaussian-process model of the high-accuracy oceanic spot data of Hoggard et al. (2017), from the coefficients they
# publish (github.com/valentineap/DynamicTopographyGP, paper_supplementary_data/sph_coeff_high_accuracy_spot.dat, md5
# f27d803b3fb627606525fd96f1fd1e0f): sum over m of (mean^2 + posterior variance). Mean square sum / 4 pi: rms 0.58 km
# water-loaded, 0.41 km air-loaded; 0.29 of the power at l 1-3, 0.39 at 4-8, 0.33 at 9-30; P30 / P2 0.035. (The
# paper's "~85% in degrees 1-3" is the share of degrees 1-15: its script sums the first 15 per-degree powers; journal 80.)
VD2020_POWER_KM2 = (0.2091, 0.6188, 0.3940, 0.3111, 0.4783, 0.2652, 0.3287, 0.2576, 0.1478, 0.1486,
                    0.1389, 0.1106, 0.0881, 0.0866, 0.0725, 0.0653, 0.0615, 0.0543, 0.0501, 0.0488,
                    0.0427, 0.0408, 0.0382, 0.0343, 0.0321, 0.0289, 0.0275, 0.0257, 0.0245, 0.0218)
VD2020_LOAD_KG_M3 = (2300., 3300.)      # the measurement's water-load and air-load density contrasts (Hoggard et al.)

# Degree correlation r_l, l = 1 .. 30, of Earth's measured residual topography (VD2020's posterior-mean coefficients,
# data/mantle/vd2020_sph_coeff_high_accuracy_spot.dat, md5 f27d803b3fb627606525fd96f1fd1e0f) with E3: Earth's own
# subduction history (Muller et al. 2019, Tectonics 38, 1884-1907) through D's unchanged SlabRecord and kernels
# (journal 70; present-day load layers of precheck-p50-geoid/e3_loads.npz with slab_alpha 'thermo', rms 645 m).
# Pass 101 (journal 136; outputs/planet-calibration/precheck-p101-degree/earth_degree_corr.py): the file's real
# harmonics read as cosine on m < 0, sine on m > 0, with the Condon-Shortley phase, the only reading that agrees with
# Holdt et al. (2022)'s independent degree 1-30 grid on the oceanic nodes (r 0.64 against <= 0.22 for every other).
# Bands: l 1-3 +0.77, 4-10 +0.07, 11-30 +0.00 (the tomographic slab proxy E1: +0.76 / +0.21 / +0.09).
EARTH_SLAB_CORR = (+0.997, +0.882, +0.322, -0.012, +0.314, +0.006, +0.030, +0.151, -0.067, +0.401,
                   -0.076, -0.107, +0.194, +0.247, +0.020, -0.045, -0.042, -0.053, +0.070, +0.175,
                   -0.101, +0.080, -0.110, -0.131, -0.076, +0.179, -0.094, -0.084, +0.028, -0.077)


def placed_share(lmax):
    """flow_placed 'earth_corr': the share of each degree's measured power the slab flow carries, max(r_l, 0)^2 from
    EARTH_SLAB_CORR (index = degree; [0] = 0; degrees beyond 30, which the measurement does not reach, 0)."""
    out = np.zeros(lmax+1)
    r = np.clip(np.asarray(EARTH_SLAB_CORR), 0., None)**2
    n = min(r.size, lmax)
    out[1:n+1] = r[:n]
    return out


def measured_power_air_m2():
    """VD2020_POWER_KM2 as an air-loaded power array indexed by degree (m^2; [0] = 0), for fields.spectrum_field."""
    w, a = VD2020_LOAD_KG_M3
    return np.concatenate(([0.], np.asarray(VD2020_POWER_KM2)*1e6*(w/a)**2))


PARAMS = {
    'dyn_std_m': 400., 'dyn_lmax': 12, 'dyn_slope': 2.5, 'dt_knot_myr': 40.,
    'hotspots_initial': 10, 'hotspot_birth_per_myr': .05,
    'hotspot_life_myr': (60., 180.),
    'hotspot_flux_km3_myr': (2000., 15000.),      # steady volcanic output (Hawaii ~1e6 km3 over ~70 Myr)
    'hotspot_radius_km': 60., 'hotspot_swell_m': 800., 'hotspot_swell_radius_km': 500.,
    'lip_volume_km3': (5e5, 4e6), 'lip_radius_km': 600., 'lip_duration_myr': 2.,
    'dyn_rule': 'random', 'rho_mantle': 3300., 'rho_water': 1030., 'threads': None,
    'slab_dip_deg': 58., 'lm_sink_mm_yr': 12., 'flow_z_cut_km': 300., 'flow_lmax': 31, 'viscosity_file': None,
    'slab_alpha': 'thermo', 'flow_power': 'flow', 'flow_placed': 'all',
}


class Mantle:
    def __init__(self, graph, rng, t0, t_end=0., params=None, field_rng=None):
        self.g, self.rng = graph, rng
        self.p = {**PARAMS, **(params or {})}
        p = self.p
        self.t0 = float(t0)
        self.knots_t = np.arange(t0-p['dt_knot_myr'], t_end+2*p['dt_knot_myr'], p['dt_knot_myr'])
        self.knots = [p['dyn_std_m']*fields.random_field(graph.xyz, p['dyn_lmax'], p['dyn_slope'], rng)
                      for _ in self.knots_t]
        self.measured = None
        self.unplaced = None
        fpl = p.get('flow_placed', 'all')
        if fpl not in ('all', 'earth_corr'):
            raise ValueError(f"flow_placed {fpl!r}: 'all' or 'earth_corr'")
        if fpl == 'earth_corr' and not (p['dyn_rule'] == 'slabs' and p.get('flow_power', 'flow') == 'measured'):
            fpl = 'all'                                     # moot without the measured slab field (default since pass 102)
        if p['dyn_rule'] == 'slabs':
            from planet_tectonics.mantle_flow import MantleFlow
            self.flow = MantleFlow(graph, p)
            if p.get('flow_power', 'flow') == 'measured':   # the start's mantle: 'measured''s field, handed over
                if field_rng is None:
                    raise ValueError("flow_power 'measured' needs its own field_rng (the start's measured field)")
                pw = measured_power_air_m2()
                self.measured = [fields.spectrum_field(graph.xyz, pw, field_rng) for _ in self.knots_t]
                if fpl == 'earth_corr':                     # the power the slabs do not account for: plate-blind,
                    sh = placed_share(pw.size-1)            # drawn after the start's field (the start stays the control's)
                    self.unplaced = [fields.spectrum_field(graph.xyz, pw*(1.-sh), field_rng) for _ in self.knots_t]
        elif p['dyn_rule'] in ('random', 'measured'):
            self.flow = None
            if p['dyn_rule'] == 'measured':
                if field_rng is None:
                    raise ValueError("dyn_rule 'measured' needs its own field_rng (the hotspots keep rng)")
                pw = measured_power_air_m2()
                self.measured = [fields.spectrum_field(graph.xyz, pw, field_rng) for _ in self.knots_t]
        else:
            raise ValueError(f"dyn_rule {p['dyn_rule']!r}: 'random', 'measured' or 'slabs'")
        # the field is held air-loaded and water-loaded under the sea (sim passes wet) with the slab flow and the
        # measured field
        self.wet_loaded = p['dyn_rule'] in ('slabs', 'measured')
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
        """Dynamic topography (m) at time t: blended random fields (dyn_rule 'random'), or the slab flow's air-loaded
        field (dyn_rule 'slabs') or the measured-spectrum fields blended at constant variance (dyn_rule 'measured'),
        x rho_m / (rho_m - rho_w) on the nodes in wet; plus hotspot swells."""
        if self.wet_loaded:
            p = self.p
            fw = 1. if wet is None else np.where(wet, p['rho_mantle']/(p['rho_mantle']-p['rho_water']), 1.)
            return self.air(t)*fw+self.swells(t)
        k = np.clip(np.searchsorted(self.knots_t, t)-1, 0, len(self.knots_t)-2)
        f = np.clip((t-self.knots_t[k])/(self.knots_t[k+1]-self.knots_t[k]), 0, 1)
        f = f*f*(3-2*f)
        out = (1-f)*self.knots[k]+f*self.knots[k+1]
        for h in self.hotspots:
            if h['start'] <= t < h['end']:
                out = out+self.p['hotspot_swell_m']*self._disc(h['pos'], self.p['hotspot_swell_radius_km']*1e3)
        return out

    def air(self, t):
        """The air-loaded field without swells (dyn_rule 'slabs' or 'measured'). With flow_power 'measured' the start's
        measured field hands over to the placed flow at constant variance over the first knot interval; with
        flow_placed 'earth_corr' the placed flow carries Earth's slab share of each degree and the plate-blind rest
        (independent of the start's field) is added to it."""
        if self.flow is not None:
            if self.measured is None:
                return self.flow.h_air
            placed = self.flow.h_air if self.unplaced is None else self.flow.h_air+self._knots_at(self.unplaced, t)
            f = float(np.clip((t-self.t0)/self.p['dt_knot_myr'], 0, 1))
            if f >= 1.:
                return placed
            a = .5*np.pi*f*f*(3-2*f)
            return np.cos(a)*self._measured_at(t)+np.sin(a)*placed
        return self._measured_at(t)

    def _measured_at(self, t):
        return self._knots_at(self.measured, t)

    def _knots_at(self, knots, t):
        k = np.clip(np.searchsorted(self.knots_t, t)-1, 0, len(self.knots_t)-2)
        f = np.clip((t-self.knots_t[k])/(self.knots_t[k+1]-self.knots_t[k]), 0, 1)
        a = .5*np.pi*f*f*(3-2*f)                          # constant variance: cos^2 + sin^2 = 1, smooth at the knots
        return np.cos(a)*knots[k]+np.sin(a)*knots[k+1]

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
