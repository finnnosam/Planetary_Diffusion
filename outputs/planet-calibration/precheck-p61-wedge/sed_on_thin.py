# Pass 61: where sediment sits on the label-free continent (column >= 20 km), by crystalline thickness under it
# (column - sediment); Earth CRUST1.0 (crsthk, sedthk) against combo-02 at 0 Myr (six seeds) and its start.
import numpy as np, glob
from planet_grid.sphere_graph import CubeSphere
from planet_earth.crust1 import C1
from planet_reviews.crust_review import column_thickness_m, run_params, sediment_bulk_m
g = CubeSphere(217); a = g.area
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
S = (np.loadtxt(C1/'sedthk.xyz')[:, 2].reshape(180, 360)*1e3)[r, c]
T = np.load('/home/claude/pre/diag/earth_thin_nodes.npz')['T'].astype(float)
B = [(0, 10), (10, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 1e9)]
def row(T, S):
    cont = T >= 20e3; X = (T-S)/1e3
    out = []
    for lo, hi in B:
        m = cont & (X >= lo) & (X < hi)
        out.append((a[m].sum()/1e12, (S[m]*a[m]).sum()/1e15))
    return np.array(out)
print('crystalline km under the label-free continent: ' + '  '.join('%g-%g' % (lo, min(hi, 99)) for lo, hi in B))
E = row(T, S)
print('Earth     area Mkm2 ' + ' '.join('%6.1f' % x for x in E[:, 0]) + '   total %.1f' % E[:, 0].sum())
print('          sed Mkm3  ' + ' '.join('%6.1f' % x for x in E[:, 1]) + '   total %.1f' % E[:, 1].sum())
for t in ('-0250.0', '+0000.0'):
    R = []
    for d in sorted(glob.glob('/home/claude/pre/c02/seed*')):
        z = np.load(f'{d}/snapshots/t{t}.npz'); p = run_params(d)
        Tp = column_thickness_m(z, p); Sp = sediment_bulk_m(z)
        R.append(row(Tp, Sp))
    R = np.mean(R, 0)
    print(f'c02 {t} area ' + ' '.join('%6.1f' % x for x in R[:, 0]) + '   total %.1f' % R[:, 0].sum())
    print(f'          sed Mkm3  ' + ' '.join('%6.1f' % x for x in R[:, 1]) + '   total %.1f' % R[:, 1].sum())
