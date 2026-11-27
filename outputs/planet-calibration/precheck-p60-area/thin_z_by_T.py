# Pass 60: elevation and area of thin crust by thickness bin, Earth (CRUST1.0 + ETOPO, ice out) and combo-02 (0 Myr, six seeds)
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_earth.reference import ice_mask
from planet_reviews.crust_review import column_thickness_m, run_params
g = CubeSphere(217); a = g.area/1e12; ice = ice_mask(g)
e = np.load('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz'); T, z = e['T'], e['z']
B = [(20, 22.5), (22.5, 25), (25, 27.5), (27.5, 30), (30, 32.5), (32.5, 35), (35, 40)]
def wmed(v, w):
    o = np.argsort(v); c = np.cumsum(w[o]); return v[o][np.searchsorted(c, c[-1]/2)] if len(v) else np.nan
print('%-12s %22s %22s %22s' % ('T bin km', 'Earth (no ice) area/zp50', 'combo-02 area/zp50', 'combo-02 z+sea-datum'))
P = []
S = 'outputs/planet-calibration/sweep-20261101-191306/combo-02'
for seed in range(1001, 1007):
    d = f'{S}/seed{seed}'; s = np.load(f'{d}/snapshots/t+0000.0.npz')
    P.append((column_thickness_m(s, run_params(d)), s['elevation_m'].astype(float), float(s['sea_level_m'])))
for lo, hi in B:
    m = ~ice & (T >= lo*1e3) & (T < hi*1e3)
    pa = []; pz = []; pz2 = []
    for Tp, zp, sl in P:
        k = (Tp >= lo*1e3) & (Tp < hi*1e3); pa.append(a[k].sum()); pz.append(wmed(zp[k], a[k])); pz2.append(wmed(zp[k]+sl, a[k]))
    print('%5.1f-%-6.1f %12.1f %9.0f %12.1f %9.0f %12.0f' % (lo, hi, a[m].sum(), wmed(z[m], a[m]), np.mean(pa), np.mean(pz), np.mean(pz2)))
print('sea level on the datum, combo-02 seeds:', [round(p[2]) for p in P])
