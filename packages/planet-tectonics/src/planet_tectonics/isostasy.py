"""Isostasy and sea level: ocean depth from age (GDH1), sediment compaction, Airy elevation, the sea level that holds
the water inventory."""
import numpy as np

from planet_tectonics.params import PARAMS



def gdh1_depth(age_myr):
    """Ocean depth (m, positive) from lithosphere age (Stein & Stein 1992)."""
    a = np.clip(age_myr, 0, None)
    return np.where(a < 20, 2600.+365.*np.sqrt(a), 5651.-2473.*np.exp(-.0278*a))



# ocean_depth_rule 'earth' / 'earth_all' (eighty-second pass, journal 107): Earth's ocean floor by age, measured on the
# n = 217 graph. ETOPO 2022 node means on dated floor (Seton et al. 2020, node >= 99% dated; CRUST1.0 total < 20 km),
# minus this module's own sediment term on GlobSed v3 (the ocean's compacting column under water) and crust term on
# CRUST1.0 (crsthk - sedthk against ocean_thick_m): the basement depth the model's floor needs to stand where Earth's
# stands under Earth's sediment and crust. 'earth': normal oceanic crust (CRUST1.0 types A1, A0), the law of cooling
# lithosphere; 'earth_all': every dated floor (plateaus and melt-affected crust in the mean: the volume of floor the
# planet has no process for, spread over all of it; read only). Area means in EARTH_FLOOR_BINS_MYR at the bin's mean
# age, pooled where an older bin is shallower (cooling lithosphere never rises as it ages; area-weighted
# pool-adjacent-violators, as lithosphere.LAMB_L_KM). GDH1 (Stein & Stein 1992: North Pacific and north-west Atlantic)
# stands 120-270 m below Earth's normal floor beyond 10 Myr (and 87 m above it on the youngest bin): normal floor
# flattens at ~5.37 km, GDH1 at 5.65. Recomputed by tests/test_ocean_depth.py (earth_floor_table; it depends on the
# ocean's compaction law and the densities)
EARTH_FLOOR_BINS_MYR = (0., 2., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 90., 100.,
                        110., 120., 130., 140., 150., 165., 200.)
EARTH_FLOOR_AGE_MYR = {
    'earth': (1.09, 3.48, 7.49, 12.48, 17.48, 22.44, 27.47, 32.47, 37.45, 42.5, 47.35, 52.51, 57.54, 62.49, 67.52, 72.46,
              77.48, 84.73, 94.96, 105.16, 114.86, 124.99, 134.64, 144.69, 157.01, 174.53),
    'earth_all': (1.1, 3.48, 7.48, 12.49, 17.47, 22.43, 27.47, 32.48, 37.45, 42.5, 47.37, 52.51, 57.54, 62.48, 67.54,
                  72.44, 77.48, 84.73, 94.98, 105.1, 114.89, 125.01, 134.77, 144.73, 156.83, 174.55)}
EARTH_FLOOR_DEPTH_M = {
    'earth': (3056.2, 3254.9, 3547.3, 3768.3, 4000.6, 4175.5, 4314.6, 4457.9, 4576.6, 4696.1, 4821.3, 4950.2, 5005.7,
              5023.2, 5071.4, 5142.6, 5150.6, 5208.7, 5346.7, 5350.1, 5350.1, 5350.1, 5413.0, 5413.0, 5413.0, 5413.0),
    'earth_all': (3059.2, 3241.7, 3534.9, 3753.1, 3991.5, 4169.6, 4297.8, 4432.8, 4547.5, 4656.2, 4741.0, 4840.0, 4903.7,
                  4903.7, 4988.9, 5068.9, 5108.5, 5147.3, 5254.8, 5254.8, 5254.8, 5254.8, 5268.1, 5268.1, 5268.1, 5272.4)}



def ocean_depth(age_myr, p=PARAMS):
    """Ocean basement depth (m, positive) from floor age under the run's ocean_depth_rule: 'gdh1' (gdh1_depth), or
    Earth's measured floor ('earth', 'earth_all': linear in sqrt(age) between the knots, flat beyond the last, and
    GDH1's own shape below the first knot, shifted to pass through it)."""
    rule = p.get('ocean_depth_rule', 'gdh1')
    if rule == 'gdh1':
        return gdh1_depth(age_myr)
    if rule not in EARTH_FLOOR_DEPTH_M:
        raise ValueError(f"ocean_depth_rule {rule!r}: 'gdh1', 'earth' or 'earth_all'")
    a = np.clip(age_myr, 0, None)
    T = np.asarray(EARTH_FLOOR_AGE_MYR[rule]); D = np.asarray(EARTH_FLOOR_DEPTH_M[rule])
    return np.where(a < T[0], D[0]-(gdh1_depth(T[0])-gdh1_depth(a)), np.interp(np.sqrt(a), np.sqrt(T), D))



def _pool_increasing(v, w):
    """Area-weighted pool-adjacent-violators: the non-decreasing sequence nearest v (weights w)."""
    blocks = [[float(x), float(y), 1] for x, y in zip(v, w)]
    i = 0
    while i < len(blocks)-1:
        if blocks[i][0] > blocks[i+1][0]:
            a, b = blocks[i], blocks.pop(i+1)
            blocks[i] = [(a[0]*a[1]+b[0]*b[1])/(a[1]+b[1]), a[1]+b[1], a[2]+b[2]]
            i = max(i-1, 0)
        else:
            i += 1
    return np.concatenate([np.full(k, x) for x, _, k in blocks])



def earth_floor_table(rule='earth', n=217, p=PARAMS, data=None, outputs=None):
    """(ages Myr, depths m) of EARTH_FLOOR_AGE_MYR / EARTH_FLOOR_DEPTH_M recomputed from the Earth references on the
    n graph (outputs/planet-reference: earth-n217.npz, earth-seafloor-age-n217.npz, globsed-n217.npz) and CRUST1.0."""
    from planet_grid.paths import DATA, OUTPUTS
    from planet_grid.sphere_graph import CubeSphere
    c1 = (data or DATA)/'crust1'; ref = (outputs or OUTPUTS)/'planet-reference'
    g = CubeSphere(n)
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(c1/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
    S = np.loadtxt(c1/'sedthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
    tn = np.array([l.split() for l in (c1/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
    z = np.load(ref/f'earth-n{n}.npz')['elevation_m'].astype(np.float64)
    ag = np.load(ref/f'earth-seafloor-age-n{n}.npz'); age = ag['age_myr'].astype(np.float64)
    gs = np.load(ref/f'globsed-n{n}.npz'); h = gs['sediment_m'].astype(np.float64)
    rc, rm, rw = p['rho_crust'], p['rho_mantle'], p['rho_water']
    floor = (T < 20e3) & (ag['dated_share'] >= .99) & (gs['valid_share'] >= .5) & (z < 0) & np.isfinite(h)
    if rule == 'earth':
        floor &= np.isin(tn, ('A1', 'A0'))
    elif rule != 'earth_all':
        raise ValueError(f"rule {rule!r}: 'earth' or 'earth_all'")
    sh, _, sm_wet = sediment_column(_solid_thickness(h, p['phi0_ocean'], p['compaction_c_ocean_per_km']),
                                    np.zeros(h.size, bool), {**p, 'sediment_compaction': True})
    depth = -(z-(sh-sm_wet/(rm-rw))-(T-S-p['ocean_thick_m'])*(rm-rc)/(rm-rw))
    a = g.area; ts, ds, ws = [], [], []
    for lo, hi in zip(EARTH_FLOOR_BINS_MYR[:-1], EARTH_FLOOR_BINS_MYR[1:]):
        m = floor & (age >= lo) & (age < hi)
        if m.sum() < 30:
            continue
        ts.append(float((age[m]*a[m]).sum()/a[m].sum())); ds.append(float((depth[m]*a[m]).sum()/a[m].sum()))
        ws.append(float(a[m].sum()))
    return np.array(ts), _pool_increasing(ds, ws)



def _solid_thickness(h, phi0, c_per_km):
    """Solid thickness S (m) of a compacting column of bulk height h (the inverse of bulk_thickness)."""
    c = c_per_km/1e3
    return h-phi0/c*(1-np.exp(-c*h))



def bulk_thickness(S, phi0, c_per_km):
    """Bulk thickness h (m) of a sediment column holding solid thickness S (m) when its porosity falls with depth z
    below the top as phi0 exp(-c z) (Athy 1930): S = h - phi0 / c (1 - exp(-c h)). Newton's method from the thin-cover
    bound S / (1 - phi0), from which it converges monotonically (the residual is convex and increasing in h)."""
    S = np.maximum(np.asarray(S, np.float64), 0.)
    phi0 = np.asarray(phi0, np.float64); c = np.asarray(c_per_km, np.float64)/1e3
    h = S/(1-phi0)
    for _ in range(50):
        e = np.exp(-c*h)
        f = h-phi0/c*(1-e)-S
        h = h-f/(1-phi0*e)
        if np.all(np.abs(f) < 1e-6):
            break
    return h



def sediment_column(sed, cont, p=PARAMS):
    """Sediment column of every node: (surface height, mass per area dry, mass per area submerged minus the water it
    displaces). With sediment_compaction, sed is solid thickness: bulk height from bulk_thickness (shale law on oceanic
    crust, shaly sandstone on continental crust), mass rho_grain S + rho_water (h - S) with water-filled pores, and
    (rho_grain - rho_water) S under water. Without it, a layer of rho_sediment: (sed, rho_sediment sed,
    (rho_sediment - rho_water) sed)."""
    rw = p['rho_water']
    if not p.get('sediment_compaction'):
        rs = p['rho_sediment']
        return sed, rs*sed, (rs-rw)*sed
    phi0 = np.where(cont, p['phi0_continent'], p['phi0_ocean'])
    c = np.where(cont, p['compaction_c_continent_per_km'], p['compaction_c_ocean_per_km'])
    h = bulk_thickness(sed, phi0, c)
    rg = p['rho_grain']
    return h, rg*sed+rw*(h-sed), (rg-rw)*sed



def isostatic_elevation(cont, thick, ocean_age, p=PARAMS):
    """Elevation (m): oceanic crust from age; continental crust from Airy thickness, water-loaded below 0 m."""
    dry = p['cont_ref_elev_m']+(thick-p['cont_ref_thick_m'])*(p['rho_mantle']-p['rho_crust'])/p['rho_mantle']
    wet = np.where(dry < 0, dry*p['rho_mantle']/(p['rho_mantle']-p['rho_water']), dry)
    return np.where(cont, wet, -ocean_depth(ocean_age, p))



def sea_level_for_volume(z, area, volume_m3):
    """Level s (m) at which the water below it, sum of (s - z) * area over z < s, equals volume_m3."""
    order = np.argsort(z)
    zs, a = z[order], area[order]
    ca = np.cumsum(a); cza = np.cumsum(zs*a)
    held = zs*ca-cza                                  # volume below level zs[k] (the k+1 lowest nodes flooded)
    k = int(np.searchsorted(held, volume_m3))
    if k == 0:
        return float(zs[0])
    k = min(k, len(zs))
    return float((volume_m3+cza[k-1])/ca[k-1])
