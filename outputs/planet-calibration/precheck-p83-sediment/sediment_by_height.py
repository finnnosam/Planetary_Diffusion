"""Pass 83 (journal 108): the continents' sediment by elevation (relative to the sea) and by the share of the
continent's area, Earth (CRUST1.0 sediment on total >= 20 km, ETOPO) against a sweep arm's 0 Myr snapshots.
Usage (repo root): python sediment_by_height.py <arm> [<arm> ...]"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
g = CubeSphere(217); a = g.area; AREA = 4*np.pi*6371e3**2
BANDS = [(-1e9, -1000), (-1000, -200), (-200, 0), (0, 250), (250, 500), (500, 1000), (1000, 2000), (2000, 3000), (3000, 1e9)]

def earth():
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
    S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    return z, T >= 20e3, S

def planet(arm, seed):
    d = np.load(f'{arm}/seed{seed}/snapshots/t+0000.0.npz')
    s = d['sediment_bulk_m'].astype(float)
    return d['elevation_m'].astype(float), (d['crust_thickness_m'].astype(float)+s) >= 20e3, s

def rows(z, cont, s):
    out = []
    for lo, hi in BANDS:
        m = cont & (z >= lo) & (z < hi)
        A = a[m].sum(); V = (s[m]*a[m]).sum()
        out.append((A/a.sum(), V*AREA/a.sum()/1e15, V/A if A else np.nan))
    return np.array(out)

if __name__ == '__main__':
    cols = [('Earth', rows(*earth()))]
    for arm in sys.argv[1:]:
        cols.append((Path(arm).name, np.mean([rows(*planet(arm, s)) for s in range(1001, 1007)], axis=0)))
    print('continent (total >= 20 km) by elevation above the sea: share of the surface / sediment Mkm3 (bulk) / mean m')
    print('band (m)        | ' + ' | '.join(f'{n:28s}' for n, _ in cols))
    for i, (lo, hi) in enumerate(BANDS):
        lab = f'{lo:>6.0f} .. {hi:<6.0f}'.replace('-1000000000', '  -inf').replace('1000000000', 'inf')
        print(f'{lab:15s} | ' + ' | '.join(f'{c[i,0]:.3f} / {c[i,1]:6.1f} / {c[i,2]:5.0f}   ' for _, c in cols))
