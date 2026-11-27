"""Pass 91: where the thin continental crust (total 20-30 km, label-free) lies, Earth against planet arms: by distance
from the continent's edge (total >= 20 km), above / under the sea, and its share of the continent in each band.
Earth: CRUST1.0 total, non-oceanic type, ice-free, ETOPO z. Args: snapshot tag, arms. Run from the repo root."""
import sys, numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import lithosphere as lt
from planet_earth.reference import ice_mask
g = CubeSphere(217)
B = [0, 50, 100, 200, 300, 500, 800, 20000]
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]
tn = np.array([l.split() for l in (DATA/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
ze = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
ice = ice_mask(g)
def table(name, tot, cont, z, keep):
    whole = cont & (tot >= 20)
    d = gl.distance_to(~whole, g.nbr, g.dist)/1e3
    a = g.area; m = whole & keep; A = a[m].sum()
    thin = m & (tot < 30)
    print(f'\n{name}: thin (20-30 km) share of continent {a[thin].sum()/A:.3f}; under the sea {a[thin & (z < 0)].sum()/max(a[thin].sum(),1):.2f}')
    print('  dist km   | cont area | thin share in band | thin area of all thin | of it under sea | thin mean z | band mean total km')
    for lo, hi in zip(B[:-1], B[1:]):
        k = m & (d >= lo) & (d < hi); kt = k & thin
        print(f'  {lo:5d}-{hi:5d} | {a[k].sum()/A:.3f} | {a[kt].sum()/max(a[k].sum(),1):.3f} | {a[kt].sum()/max(a[thin].sum(),1):.3f} | '
              f'{a[kt & (z<0)].sum()/max(a[kt].sum(),1):.2f} | {np.average(z[kt], weights=a[kt]) if kt.any() else np.nan:7.0f} | '
              f'{np.average(tot[k], weights=a[k]) if k.any() else np.nan:5.1f}')
cont_e = ~np.isin(tn, lt.OCEANIC_TYPES)
table('Earth (CRUST1.0, ice-free)', T, cont_e, ze, ~ice)
tag = sys.argv[1]
def summary(arm, tag):
    """Six-seed means per band: continent area share, thin share in band, band mean total."""
    rows = []
    for s in range(1001, 1007):
        q = np.load(f'{arm}/seed{s}/snapshots/{tag}.npz')
        tot = (q['crust_thickness_m']+q['sediment_m'])/1e3; cont = q['continental'].astype(bool)
        whole = cont & (tot >= 20); d = gl.distance_to(~whole, g.nbr, g.dist)/1e3; a = g.area; A = a[whole].sum()
        thin = whole & (tot < 30); r = [a[thin].sum()/A]
        for lo, hi in zip(B[:-1], B[1:]):
            k = whole & (d >= lo) & (d < hi)
            r += [a[k].sum()/A, a[k & thin].sum()/max(a[k].sum(), 1), np.average(tot[k], weights=a[k]) if k.any() else np.nan]
        rows.append(r)
    m = np.nanmean(rows, 0)
    print(f'\n{arm.split("/")[-1]} {tag} (six-seed mean): thin share {m[0]:.3f}')
    for i, (lo, hi) in enumerate(zip(B[:-1], B[1:])):
        print(f'  {lo:5d}-{hi:5d} | {m[1+3*i]:.3f} | {m[2+3*i]:.3f} | {m[3+3*i]:5.1f}')
for arm in sys.argv[2:]:
    if '@' in arm: arm, t = arm.split('@')
    else: t = tag
    summary(arm, t)
