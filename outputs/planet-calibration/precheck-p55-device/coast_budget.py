"""Pass 55: land by distance inland (coast_profile's bins), per band: z on the datum, crystalline crust, sediment, the
planet's Airy line (cont_ref 38 km at +300 m, rho_c 2,830, rho_m 3,300; sediment at 2,400 and, for the planet, at its
own bulk density from compaction ~ printed beside), and z - Airy. Distance bins do not depend on CRUST1.0's thickness
errors, so band medians of z - Airy are free of the regression dilution that biases thickness bins (pass 34).
The change of z - Airy from the coast band to the interior band is what the lithosphere (C) and anything else not in
the planet's Airy line does to the rise."""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import DATA, OUTPUTS
from planet_reviews.crust_review import sediment_bulk_m
from planet_reviews.coast_profile import distance_to_ocean_km, BINS_KM
from planet_earth.reference import ice_mask
RC, RM, RS = 2830., 3300., 2400.
def wmed(v, w):
    if w.sum() == 0: return np.nan
    o = np.argsort(v); c = np.cumsum(w[o]); return float(v[o][np.searchsorted(c, c[-1]/2)])
def airy(Tc, S): return 300. + (Tc-38000.)*(1-RC/RM) + S*(1-RS/RM)
g = CubeSphere(217); a = g.area
lat = np.degrees(np.arcsin(g.xyz[:, 2])); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
ld = lambda f: np.loadtxt(DATA/'crust1'/f)[:, 2].reshape(180, 360)[ii, jj]*1e3
T, S = ld('crsthk.xyz'), ld('sedthk.xyz')
ze = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
keys = ['z', 'Tc', 'S', 'airy', 'resid']
def bands(zd, Tc, S, land, d):
    out = []
    for lo, hi in zip(BINS_KM[:-1], BINS_KM[1:]):
        m = land & (d >= lo) & (d < hi); w = a[m]
        out.append([wmed(zd[m], w), wmed(Tc[m], w)/1e3, wmed(S[m], w), wmed(airy(Tc, S)[m], w), wmed((zd-airy(Tc, S))[m], w)])
    return np.array(out)
def show(name, b):
    print(name)
    for (lo, hi), r in zip(zip(BINS_KM[:-1], BINS_KM[1:]), b):
        print(f'   {lo:>5}-{hi:<5} km: z {r[0]:6.0f}  crystalline {r[1]:5.1f} km  sediment {r[2]:5.0f} m  Airy {r[3]:6.0f}  z - Airy {r[4]:6.0f}')
    print(f'   rise (600-1200 minus 0-50): z {b[4,0]-b[0,0]:.0f} m, Airy {b[4,3]-b[0,3]:.0f} m, z - Airy {b[4,4]-b[0,4]:.0f} m')
landE = (ze > 0) & ~ice_mask(g)
show('Earth (sea level = datum)', bands(ze, T-S, S, landE, distance_to_ocean_km(g, landE)))
for combo in sys.argv[2:]:
    bs = []
    for sd in sorted((Path(sys.argv[1])/combo).glob('seed*')):
        z = np.load(sd/'snapshots'/'t+0000.0.npz'); e = z['elevation_m'].astype(float)
        Sb = np.where(z['continental'], sediment_bulk_m(z), 0.)
        bs.append(bands(e+float(z['sea_level_m']), z['crust_thickness_m'].astype(float), Sb, e > 0, distance_to_ocean_km(g, e > 0)))
    show(f'{combo} (6-seed mean of band medians, z on the datum)', np.nanmean(bs, 0))
