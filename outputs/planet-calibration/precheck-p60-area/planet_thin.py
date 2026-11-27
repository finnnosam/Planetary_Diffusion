# Pass 60: where the planet's thin continental crust (crust + sediment 20-30 km, label-free) lies, combo-02 and the
# control (combo-04) at 0 Myr, six seeds; same cuts as earth_thin.py.
import json, sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, run_params
S = 'outputs/planet-calibration/sweep-20261101-191306'
g = CubeSphere(217); a = g.area/1e12; nb = g.nbr; ok = nb >= 0
def cuts(f, params):
    z = np.load(f); T = column_thickness_m(z, params); flag = z['continental'].astype(bool)
    e = z['elevation_m'].astype(float)                       # relative to sea level already
    cont = T >= 20e3; thin = cont & (T < 30e3); bnd = z['boundary_type'] != 0
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    de = gl.distance_to(edge, nb, g.dist, 3000e3)/1e3; db = gl.distance_to(bnd, nb, g.dist, 5000e3)/1e3
    la = z['lithosphere_age_myr'].astype(float); ca = z['crust_age_myr'].astype(float)
    o = {'cont': a[cont].sum(), 'thin': a[thin].sum(), 'thin_flag': a[thin & flag].sum(), 'thin_oceanflag': a[thin & ~flag].sum(),
         'thin_sub': a[thin & (e < 0)].sum(), 'cont_sub': a[cont & (e < 0)].sum(),
         'cont_flag_lt20': a[flag & ~cont].sum(),
         'thin_lt150_edge': a[thin & (de < 150)].sum(), 'thin_150_500': a[thin & (de >= 150) & (de < 500)].sum(), 'thin_gt500': a[thin & (de >= 500)].sum(),
         'thin_lt300_bnd': a[thin & (db < 300)].sum(),
         'thin_lith_lt100': a[thin & (la < 100)].sum(), 'thin_litharc': np.nan,
         'edge_nodes': int(edge.sum())}
    for lo, hi in ((-1e9, -2000), (-2000, -200), (-200, 0), (0, 500), (500, 1e9)):
        o[f'thin_z{int(max(lo,-9999))}'] = a[thin & (e >= lo) & (e < hi)].sum()
    return o
for combo in ('combo-02', 'combo-04'):
    R = []
    for seed in range(1001, 1007):
        d = f'{S}/{combo}/seed{seed}'
        R.append(cuts(f'{d}/snapshots/t+0000.0.npz', run_params(d)))
    print(combo, '(six-seed means, Mkm2):')
    for k in R[0]:
        print('   %-18s %8.2f   (%s)' % (k, np.mean([r[k] for r in R]), ' '.join('%.1f' % r[k] for r in R)))
