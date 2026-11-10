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

Frames (lith_frame; journal 104). 'artemieva' (above): Artemieva's thermal-age L with Lamb's central drho 30, two
frames joined (journal 102 section 4). 'lamb': L measured from a global seismic lithosphere model on Earth's own
ages, with the density excess Lamb et al. derived for seismic L. Hoggard et al. 2020 (Nat. Geosci. 13, 504;
https://osf.io/twksd), SL2013sv LAB depth, sampled at the n = 217 nodes of Earth's continent as above (earth_continent),
crust <= 45 km (unthickened: the planet thickens collided lithosphere itself, collision_lithosphere 'thicken'); the
area mean of each age class (CLASS_AGE_MA), pooled where an older class is thinner than a younger one (conductive
cooling never thins lithosphere as its thermal age grows; area-weighted pool-adjacent-violators), linear between the
classes, flat below the youngest class and beyond the oldest (LAMB_AGE_MA, LAMB_L_KM): 94 km at 33 Ma, 112 by 100
Ma, 114 through the late Proterozoic, 165 at 1.5 Ga, 178 for cratons (Artemieva's line: 95 / 98 / 124 / 154 / 164-222).
Its density excess is Lamb's k (10.753 m of elevation per km of L, the reduction of his ds02 table) x rho_mantle:
35.5 kg/m^3; its reference is the table's area mean on Earth's whole continent (130.68 km); lith_max_km is the
table's top (177.74 km). The conductive-cooling form L_inf - (L_inf - L_a) exp(-t / tau) does not describe these
means (tau runs to 1,500 Myr, the fit's bound; rms 54 km; precheck-p80-hoggard/hoggard_law.txt). The table rises by
LAMB_EPS km per Myr besides, so it inverts (stretching, thickening); 0.004 km over 3.6 Gyr.
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


# lith_frame 'lamb' (module text; recomputed by tests/test_lithosphere.py when data/Hoggard2020 is present)
LAMB_AGE_MA = (0., 33., 100., 270., 770., 1520., 1750., 3600.)
LAMB_L_KM = (93.78, 93.78, 111.74, 114.29, 114.29, 164.50, 177.74, 177.74)
LAMB_EPS = 1e-6
LAMB_REF_KM = 130.68
LAMB_K_M_PER_KM = 10.753   # Lamb et al. 2020, ds02's reduced elevation z + k (L - 100)
HOGGARD_SL2013SV = ('Hoggard2020', 'twksd-osfstorage-archive', 'supplementary_datasets', 'lithospheric_thickness_maps',
                    'global', 'SL2013sv', 'SL2013sv.csv')


def _lamb(p):
    return p.get('lith_frame', 'artemieva') == 'lamb'


def frame(p):
    """(reference L km, density excess kg/m^3, top L km) of the run's lithosphere frame."""
    if _lamb(p):
        return LAMB_REF_KM, LAMB_K_M_PER_KM*p['rho_mantle']/1e3, LAMB_L_KM[-1]+LAMB_EPS*LAMB_AGE_MA[-1]
    return p['lith_ref_km'], p['lith_drho_kg_m3'], p['lith_max_km']


def thickness_km(age_myr, p):
    if _lamb(p):
        t = np.maximum(age_myr, 0.)
        return np.interp(t, LAMB_AGE_MA, LAMB_L_KM)+LAMB_EPS*np.minimum(t, LAMB_AGE_MA[-1])
    return p['lith_L0_km']+p['lith_rate_km_myr']*np.maximum(age_myr, 0.)


def age_of_thickness(L, p):
    """The thermal age (Myr, >= 0) whose lithosphere is L km thick (the inverse of thickness_km)."""
    if _lamb(p):
        tk = np.asarray(LAMB_AGE_MA)
        return np.interp(L, np.asarray(LAMB_L_KM)+LAMB_EPS*tk, tk)
    return np.maximum((L-p['lith_L0_km'])/p['lith_rate_km_myr'], 0.)


def earth_lamb_table(n=217, data=None):
    """(LAMB_L_KM at the class ages, LAMB_REF_KM) recomputed from Hoggard et al. 2020's SL2013sv map and CRUST1.0
    (see the module text)."""
    from planet_grid.paths import DATA
    from planet_grid.sphere_graph import CubeSphere
    from planet_earth.reference import ice_mask
    root = data or DATA
    g = CubeSphere(n)
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(root/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]
    tn = np.array([l.split() for l in (root/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
    cont = ~np.isin(tn, OCEANIC_TYPES) & (T >= 20) & ~ice_mask(g)
    d = np.loadtxt(root.joinpath(*HOGGARD_SL2013SV), delimiter=',')
    grid = np.full((361, 721), np.nan)
    grid[np.round((90-d[:, 1])*2).astype(int), np.round((d[:, 0] % 360)*2).astype(int) % 720] = d[:, 2]
    grid[:, 720] = grid[:, 0]
    L = grid[np.round((90-lat)*2).astype(int), np.round((lon % 360)*2).astype(int)]
    age = np.array([CLASS_AGE_MA.get(t, np.nan) for t in tn])
    m = cont & (T <= 45) & np.isfinite(L)
    ages = np.array(sorted(set(age[m])))
    mu = [float((L[m & (age == t)]*g.area[m & (age == t)]).sum()/g.area[m & (age == t)].sum()) for t in ages]
    wa = [float(g.area[m & (age == t)].sum()) for t in ages]
    blocks = [[v, w, [i]] for i, (v, w) in enumerate(zip(mu, wa))]   # pool adjacent violators, area-weighted
    i = 0
    while i < len(blocks)-1:
        if blocks[i][0] > blocks[i+1][0]:
            a_, b_ = blocks[i], blocks[i+1]
            blocks[i] = [(a_[0]*a_[1]+b_[0]*b_[1])/(a_[1]+b_[1]), a_[1]+b_[1], a_[2]+b_[2]]; del blocks[i+1]
            i = max(i-1, 0)
        else:
            i += 1
    iso = np.empty(len(mu))
    for v, _, idx in blocks:
        iso[idx] = v
    whole = np.interp(age[cont], ages, iso)
    return dict(zip(ages.tolist(), iso.tolist())), float((whole*g.area[cont]).sum()/g.area[cont].sum())


def elevation_term(age_myr, p):
    """Surface offset (m, dry) of a continental column whose lithosphere has this tectonothermal age (NaN -> 0)."""
    L = thickness_km(np.nan_to_num(age_myr, nan=0.), p)
    ref, drho, _ = frame(p)
    dz = -drho/p['rho_mantle']*(L-ref)*1e3
    return np.where(np.isnan(age_myr), 0., dz)


def stretched_age(age_myr, ratio, p):
    """Age after uniform stretching that leaves the crust at `ratio` (T_after / T_before) of its thickness: the
    lithosphere thins by the same factor and keeps the thermal age that thickness corresponds to (>= 0)."""
    return age_of_thickness(thickness_km(age_myr, p)*ratio, p)


def thickened_therm(therm, thick, gain, t, p):
    """collision_lithosphere 'thicken' (journal 101): continental nodes whose crust gains `gain` (m) on `thick` (m)
    by collision take the thermal age whose L is L (thick + gain) / thick (pure shear, the mirror of stretched_age).
    Returns (new therm, changed mask); NaN therm (no lithosphere age) and gain <= 0 are left alone.
    Thickened lithosphere beyond the thickest stable cratonic root is convectively unstable and removed (Houseman,
    McKenzie & Molnar 1981): L is held at lith_max_km (C's own law at the oldest crust, 3.6 Ga: 237.6 km), and a node
    already thicker keeps its L. Without the bound pure shear ratchets: collision volume arrives step after step
    while the collapse flow and foundering hold the crust near 50-65 km (journal 101: L > 300 km on 1% of the
    continent by -40 Myr, then without bound)."""
    therm = np.asarray(therm, np.float64).copy()
    m = (np.asarray(gain) > 0) & np.isfinite(therm)
    if m.any():
        ratio = (thick[m]+gain[m])/np.maximum(thick[m], 1e3)
        L0 = thickness_km(t-therm[m], p)
        L1 = np.minimum(L0*ratio, np.maximum(L0, frame(p)[2]))
        therm[m] = t-age_of_thickness(L1, p)
    return therm, m
