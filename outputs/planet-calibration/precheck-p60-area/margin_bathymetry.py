# Pass 60: Earth's ocean depth seaward of passive continental edges on n = 217 (ETOPO, ice out): p25/50/75 by distance
# from the label-free continent (CRUST1.0 >= 20 km); ocean nodes > 300 km from any Bird boundary. combo-02 beside it
# (its own boundaries, label-free continent, six seeds at 0 Myr; elevation relative to its sea level).
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_earth.bird import earth_boundaries
from planet_earth.reference import ice_mask
from planet_reviews.crust_review import column_thickness_m, run_params
g = CubeSphere(217); a = g.area; nb = g.nbr
DB = [(0, 50), (50, 100), (100, 150), (150, 200), (200, 300), (300, 500), (500, 700), (700, 1000), (1000, 1500)]
def wq(v, w):
    o = np.argsort(v); c = np.cumsum(w[o]); c /= c[-1]; return [v[o][min(np.searchsorted(c, q), len(v)-1)] for q in (.25, .5, .75)]
def prof(z, cont, bnd, keep):
    dc = gl.distance_to(cont, nb, g.dist, 3e6)/1e3; far = gl.distance_to(bnd, nb, g.dist, 5e6)/1e3 > 300
    out = []
    for lo, hi in DB:
        m = ~cont & far & keep & (dc > lo) & (dc <= hi)
        out.append(wq(z[m], a[m]) + [a[m].sum()/1e12] if m.any() else [np.nan]*4)
    return np.array(out, float)
e = np.load('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz')
E = prof(e['z'].astype(float), e['T'] >= 20e3, earth_boundaries(g), ~ice_mask(g))
Ps = []
for seed in range(1001, 1007):
    d = f'outputs/planet-calibration/sweep-20261101-191306/combo-02/seed{seed}'; s = np.load(f'{d}/snapshots/t+0000.0.npz')
    Ps.append(prof(s['elevation_m'].astype(float), column_thickness_m(s, run_params(d)) >= 20e3, s['boundary_type'] != 0, np.ones(g.size, bool)))
Pm = np.nanmean(Ps, 0)
print('km seaward    Earth z p25/50/75 (area Mkm2)        combo-02 z p25/50/75 (area)')
for (lo, hi), r, q in zip(DB, E, Pm):
    print('%5d-%-5d  %6.0f %6.0f %6.0f (%5.1f)     %6.0f %6.0f %6.0f (%5.1f)' % (lo, hi, *r, *q))
