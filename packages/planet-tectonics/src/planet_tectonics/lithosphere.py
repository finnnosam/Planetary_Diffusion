"""C: the continental lithosphere's thermal state (fifty-sixth pass, journal 78; designed in journal 77 section 6).

Each continental node carries the time of its last thermal event (plate raster `therm`, Myr; NaN on oceanic crust).
Its tectonothermal age t (Myr) sets the thermal lithosphere's thickness by Artemieva's (2006, TC1, Tectonophysics 416,
245-277) regression over continental geotherms from > 3.6 Ga to the present,

    L = lith_L0_km + lith_rate_km_myr t = 93.6 + 0.04 t   (km),

and the mantle lithosphere, denser than the asthenosphere by lith_drho_kg_m3 (Lamb et al. 2020, G-cubed, "Global whole
lithosphere isostasy": 20-40 kg/m^3 net of depletion; 30 the central value), stands the surface at

    dz = -(drho / rho_mantle) (L - lith_ref_km)

relative to the planet's Airy line (locally compensated: the load is far broader than the flexural wavelength;
x rho_m / (rho_m - rho_w) under the sea). lith_ref_km is Earth's area mean of L over its continent from the same age
map (EARTH_LITH_REF_KM), so on Earth's own age distribution the term moves no volume: it lowers old interiors and
raises young margins about the datum the planet's Airy line was set on.

Ages from Earth (start): every CRUST1.0 crust type (CNtype1, Laske et al. 2013) takes the midpoint of the geological
time-scale interval its name gives (CLASS_AGE_MA); the planet's start continent is joined to Earth's continent by
graph distance from the continent's edge, in km (EARTH_LITH_AGE_TABLE: Earth's age quantiles in each distance bin
DIST_BINS_KM, the last pooling every node >= 1,000 km inland; Earth's continents have young rims and old cores: the
term averages +271 m within 50 km of the edge and -135 m beyond 1,000 km, as accreted belts have a width, not a
share of the continent), and within a bin by the rank of the start's own craton age, so each start craton keeps one
place in its bin's age distribution (old cratons old, coherent in space). The first builds joined by crystalline
thickness rank (little of Earth's age structure) and by distance rank (a supercontinent's rim takes all of Earth's
young bins); the pre-checks' interiors came out young (journal 78). tests/test_lithosphere.py recomputes the table
from data/crust1.

Resets (PlateSim): age 0 at a collision contact (crust consumed under continental crust, or continental nodes on a
collision boundary that is closing faster than boundary_min_mm_yr), in the continental arc belt over a converging
trench, on new continental crust from arcs; ocean-floor sediment basins that turn continental keep the floor's age.
Stretching thins the lithosphere with the crust (McKenzie 1978 uniform stretching): L -> L T_after / T_before; new
rift crust starts at age 0.
"""
import numpy as np


# tectonothermal age (Ma) of each continental CRUST1.0 type: the midpoint of the time-scale interval in its name
# (Archean 2.5-4.0 Ga with crust from 3.6; early Archean 3.6-2.8, late 2.8-2.5; early/mid Proterozoic 2.5-1.0; late
# Proterozoic 1.0-0.541; Phanerozoic 0.541-0; platforms: Precambrian basement below the Archean, 2.5-0.541; active
# orogens, arcs, rifts and extended crust: Cenozoic, 66-0; margins, shelves and transitions: rifted since Pangea's
# breakup, 200-0)
CLASS_AGE_MA = {'G1': 3200., 'G2': 2650., 'F-': 3050., 'H1': 1750., 'H2': 1750., 'I1': 770., 'I2': 770.,
                'D-': 1520., 'E-': 1520., 'Z1': 270., 'Z2': 270., 'Y1': 270.,
                'T-': 100., 'U-': 100., 'C-': 100., 'S-': 100., 'W-': 100., 'Y2': 100.,
                **{c: 33. for c in ('O-', 'P-', 'Q-', 'R1', 'R2', 'J-', 'K-', 'L1', 'L2', 'M-', 'N-', 'X-', 'V2')}}
OCEANIC_TYPES = ('A1', 'A0', 'B-', 'V1', 'Y3')
DIST_BINS_KM = (0., 50., 100., 150., 200., 300., 400., 500., 600., 800., 1000., float('inf'))
TABLE_V = tuple(k/10 for k in range(11))   # age quantiles within a bin


def earth_continent(n=217, data=None):
    """Earth's continent on the n graph (CRUST1.0 type continental, total >= 20 km, ice sheets excluded, as the
    planet_reviews Earth references): (area, graph distance to non-continental crust m (ice sheets' crust counted as
    continent), tectonothermal age Ma)."""
    from planet_grid.paths import DATA
    from planet_grid.sphere_graph import CubeSphere
    from planet_earth.reference import ice_mask
    from planet_landscape import graph_landscape as gl
    c1 = (data or DATA)/'crust1'
    g = CubeSphere(n)
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(c1/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
    tn = np.array([l.split() for l in (c1/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
    whole = ~np.isin(tn, OCEANIC_TYPES) & (T >= 20e3)
    d = gl.distance_to(~whole, g.nbr, g.dist)
    cont = whole & ~ice_mask(g)
    age = np.array([CLASS_AGE_MA[t] for t in tn[cont]])
    return g.area[cont], d[cont], age


def _wquantile(v, w, q):
    o = np.argsort(v, kind='stable'); cw = np.cumsum(w[o]); cw = cw/cw[-1]
    return v[o][np.minimum(np.searchsorted(cw, np.asarray(q)), v.size-1)]


def earth_table(n=217, data=None):
    """(EARTH_LITH_REF_KM, EARTH_LITH_AGE_TABLE) recomputed from CRUST1.0 (see the module text)."""
    a, d, age = earth_continent(n, data)
    L = 93.6+.04*age
    ref = float((L*a).sum()/a.sum())
    rows = []
    for lo, hi in zip(DIST_BINS_KM[:-1], DIST_BINS_KM[1:]):
        m = (d >= lo*1e3) & (d < hi*1e3)
        rows.append([float(x) for x in _wquantile(age[m], a[m], TABLE_V)])
    return ref, rows


# recomputed by tests/test_lithosphere.py (n = 217)
EARTH_LITH_REF_KM = 130.67
# rows: DIST_BINS_KM, distance from the continent's edge (edge first); columns: age quantiles TABLE_V (Ma)
EARTH_LITH_AGE_TABLE = (
    (33., 33., 100., 100., 100., 100., 100., 100., 100., 100., 3200.),
    (33., 33., 33., 100., 100., 100., 100., 100., 100., 770., 3200.),
    (33., 33., 33., 100., 100., 100., 100., 100., 770., 770., 3200.),
    (33., 33., 33., 33., 100., 100., 100., 270., 770., 1750., 3200.),
    (33., 33., 33., 33., 100., 100., 270., 770., 1520., 1750., 3200.),
    (33., 33., 33., 100., 100., 270., 770., 770., 1520., 2650., 3200.),
    (33., 33., 33., 100., 100., 270., 770., 770., 1750., 2650., 3200.),
    (33., 33., 33., 100., 100., 270., 770., 1520., 1750., 2650., 3200.),
    (33., 33., 33., 100., 270., 770., 1520., 1750., 1750., 3200., 3200.),
    (33., 33., 100., 770., 770., 1520., 1750., 1750., 1750., 3200., 3200.),
    (33., 33., 33., 770., 770., 1520., 1750., 1750., 1750., 2650., 3200.),
)


def start_ages(dist, v, table=None):
    """Tectonothermal ages (Ma) for start continental nodes: the distance from the continent's edge (m) picks the
    row of Earth's table (DIST_BINS_KM), the craton rank v in [0, 1] the age quantile within it (linear between the
    table's quantiles)."""
    tab = np.asarray(EARTH_LITH_AGE_TABLE if table is None else table, np.float64)
    b = np.clip(np.searchsorted(DIST_BINS_KM, np.asarray(dist)/1e3, side='right')-1, 0, len(DIST_BINS_KM)-2)
    x = np.clip(np.asarray(v, np.float64), 0., 1.)*(len(TABLE_V)-1)
    k = np.clip(x.astype(int), 0, len(TABLE_V)-2); f = x-k
    return tab[b, k]*(1-f)+tab[b, k+1]*f


def thickness_km(age_myr, p):
    return p['lith_L0_km']+p['lith_rate_km_myr']*np.maximum(age_myr, 0.)


def elevation_term(age_myr, p):
    """Surface offset (m, dry) of a continental column whose lithosphere has this tectonothermal age (NaN -> 0)."""
    L = thickness_km(np.nan_to_num(age_myr, nan=0.), p)
    dz = -p['lith_drho_kg_m3']/p['rho_mantle']*(L-p['lith_ref_km'])*1e3
    return np.where(np.isnan(age_myr), 0., dz)


def stretched_age(age_myr, ratio, p):
    """Age after uniform stretching that leaves the crust at `ratio` (T_after / T_before) of its thickness: the
    lithosphere thins by the same factor and keeps the thermal age that thickness corresponds to (>= 0)."""
    L = thickness_km(age_myr, p)*ratio
    return np.maximum((L-p['lith_L0_km'])/p['lith_rate_km_myr'], 0.)
