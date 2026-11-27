"""Pass 103 probe: how clustered the continents are, Earth against the planet: the label-free continent's indicator
function (crust + bulk sediment >= 20 km; Earth CRUST1.0) expanded to l 1-10, its power share at l 1, 2, 3 and 4-10
of l 1-10, and the resultant length of the continent (|area-weighted mean of the unit vectors| / continental area:
1 = all in one point, 0 = spread evenly). Per arm and seed at 0 Myr. Read. Usage: python3 cont_spectrum.py <arm> ..."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_reviews.crust_review import sediment_bulk_m
g = CubeSphere(217); a = g.area; LM = 10; L = mf.alm_l(LM); M = mf.alm_m(LM); w = np.where(M == 0, 1., 2.)
def row(c):
    C = mf.point_alm(g.xyz, (c-c.mean()).astype(float)*a/g.radius**2, LM)
    P = np.bincount(L, w*np.abs(C)**2, LM+1)[1:]; P = P/P.sum()
    v = (a[c, None]*g.xyz[c]).sum(0); return P[0], P[1], P[2], P[3:].sum(), np.linalg.norm(v)/a[c].sum()
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
print('                     share of l 1-10 power: l1 / l2 / l3 / l4-10 | resultant')
print('Earth (CRUST1.0)     %.2f / %.2f / %.2f / %.2f | %.2f' % row((T+S) >= 20e3))
for arm in map(Path, sys.argv[1:]):
    R = []
    for s in sorted(int(p.name[4:]) for p in arm.glob('seed*')):
        for t in ('t+0000.0', 't-0250.0'):
            q = np.load(arm/f'seed{s}'/'snapshots'/f'{t}.npz')
            R.append(row((q['crust_thickness_m']+sediment_bulk_m(q)) >= 20e3))
    R = np.array(R)
    print('%-20s %.2f / %.2f / %.2f / %.2f | %.2f   (by seed: %s)' % ((arm.name[:20],)+tuple(R[0::2].mean(0))+(', '.join('%.2f' % r[4] for r in R[0::2]),)))
    print('%-20s %.2f / %.2f / %.2f / %.2f | %.2f' % (('  the start, -250',)+tuple(R[1::2].mean(0))))
