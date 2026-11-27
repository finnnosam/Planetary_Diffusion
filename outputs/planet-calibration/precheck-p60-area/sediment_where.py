# Pass 60: where sediment lies. Earth: CRUST1.0 sedthk (and GlobSed on its valid ocean area) split by label-free
# continent (crsthk >= 20 km) / ocean; planet combo-02 and combo-04 at 0 Myr (bulk sediment, six seeds).
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_earth.crust1 import C1
from planet_reviews.crust_review import column_thickness_m, run_params, sediment_bulk_m
g = CubeSphere(217); a = g.area
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
S = (np.loadtxt(C1/'sedthk.xyz')[:, 2].reshape(180, 360)*1e3)[r, c]
T = np.load('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz')['T'].astype(float)
cont = T >= 20e3
gs = np.load('outputs/planet-reference/globsed-n217.npz'); G, vs = gs['sediment_m'], gs['valid_share']
V = lambda m, s: (s[m]*a[m]).sum()/1e15   # Mkm3
print('Earth CRUST1.0 sediment: on continent %.0f Mkm3, on ocean (< 20 km) %.0f Mkm3; ocean mean %.0f m over %.0f Mkm2'
      % (V(cont, S), V(~cont, S), V(~cont, S)/a[~cont].sum()*1e15, a[~cont].sum()/1e12))
ok = ~cont & (vs > .5) & np.isfinite(G)
print('Earth GlobSed on ocean nodes with valid data (%.0f Mkm2): %.0f Mkm3 (mean %.0f m); CRUST1.0 on the same nodes %.0f Mkm3'
      % (a[ok].sum()/1e12, V(ok, G), V(ok, G)/a[ok].sum()*1e15, V(ok, S)))
for combo in ('combo-02', 'combo-04'):
    o = []
    for seed in range(1001, 1007):
        d = f'outputs/planet-calibration/sweep-20261101-191306/{combo}/seed{seed}'; z = np.load(f'{d}/snapshots/t+0000.0.npz')
        Tp = column_thickness_m(z, run_params(d)); sb = sediment_bulk_m(z).astype(float); cp = Tp >= 20e3
        o.append((V(cp, sb), V(~cp, sb), V(~cp, sb)/a[~cp].sum()*1e15))
    o = np.array(o); print('%s sediment: on continent %.0f, on ocean %.0f Mkm3 (ocean mean %.0f m)  [seeds ocean: %s]' % (combo, *o.mean(0), ' '.join('%.0f' % x for x in o[:, 1])))
