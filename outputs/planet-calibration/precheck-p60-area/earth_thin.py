# Pass 60: where Earth's thin continental crust (CRUST1.0 total 20-30 km, label-free) lies: CRUST1.0 type, elevation
# (ETOPO on n = 217), distance from a passive edge / from a plate boundary.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.margin_width import earth_nodes, BOUNDARY_KM
from planet_earth.bird import earth_boundaries
from planet_earth.crust1 import C1
from planet_earth.reference import ice_mask
g = CubeSphere(217); a = g.area/1e12
T = np.loadtxt(C1/'crsthk.xyz')[:, 2].reshape(180, 360)*1e3
types = np.array([l.split() for l in (C1/'CNtype1-1.txt').read_text().splitlines() if l.strip()])
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
T = T[r, c]; ty = types[r, c]
z = np.load('outputs/planet-reference/earth-n217.npz')['elevation_m']
GROUPS = {'shelf C-': ['C-'], 'slope/transition S- T- U-': ['S-', 'T-', 'U-'], 'extended M- N- X- V2 Y1 Y2': ['M-', 'N-', 'X-', 'V2', 'Y1', 'Y2'],
          'arc/forearc J- K- L1 L2': ['J-', 'K-', 'L1', 'L2'], 'oceanic plateau B- W-': ['B-', 'W-'],
          'oceanic A0 A1 V1 Y3': ['A0', 'A1', 'V1', 'Y3'],
          'orogen O- P- Q- R1 R2': ['O-', 'P-', 'Q-', 'R1', 'R2'],
          'shield/platform/Phanerozoic': ['D-', 'E-', 'F-', 'G1', 'G2', 'H1', 'H2', 'I1', 'I2', 'Z1', 'Z2']}
cont = T >= 20e3; thin = cont & (T < 30e3); A = a[cont].sum()
bnd = earth_boundaries(g)
nb = g.nbr; ok = nb >= 0
edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
far = gl.distance_to(bnd, nb, g.dist, 5000e3) > BOUNDARY_KM*1e3
dpe = gl.distance_to(edge & far, nb, g.dist, 3000e3)/1e3
de = gl.distance_to(edge, nb, g.dist, 3000e3)/1e3
db = gl.distance_to(bnd, nb, g.dist, 5000e3)/1e3
ice = ice_mask(g)
print('Earth label-free continent %.1f Mkm2 (%.3f of surface); thin 20-30 km %.1f Mkm2 (%.3f of continent, %.3f of surface)'
      % (A, A/a.sum(), a[thin].sum(), a[thin].sum()/A, a[thin].sum()/a.sum()))
print('%-30s %8s %8s %8s %9s %9s %9s' % ('CRUST1.0 type group', 'cont', 'thin', 'thin/grp', 'z<0 thin', 'z p50', 'near bnd'))
for k, v in GROUPS.items():
    m = cont & np.isin(ty, v); mt = m & thin
    zz = z[mt]; o = np.argsort(zz); w = np.cumsum(a[mt][o])
    zp = zz[o][np.searchsorted(w, w[-1]/2)] if mt.any() else np.nan
    print('%-30s %8.1f %8.1f %8.2f %9.2f %9.0f %9.2f' % (k, a[m].sum(), a[mt].sum(), a[mt].sum()/max(a[m].sum(), 1e-9),
          a[mt & (z < 0)].sum()/max(a[mt].sum(), 1e-9), zp, a[mt & (db < 300)].sum()/max(a[mt].sum(), 1e-9)))
print('thin crust by elevation: z<-2000 %.1f, -2000..-200 %.1f, -200..0 %.1f, 0..500 %.1f, >500 %.1f Mkm2' % tuple(
    a[thin & (z >= lo) & (z < hi)].sum() for lo, hi in ((-1e9, -2000), (-2000, -200), (-200, 0), (0, 500), (500, 1e9))))
print('thin crust submerged: %.3f; continent submerged (z<0): %.3f; ice nodes in thin: %.2f Mkm2' % (
    a[thin & (z < 0)].sum()/a[thin].sum(), a[cont & (z < 0)].sum()/A, a[thin & ice].sum()))
for lab, dd in (('distance from any continental edge', de), ('distance from a passive edge (>300 km from a boundary)', dpe), ('distance from a plate boundary', db)):
    print(lab + ': ' + ', '.join('%d-%d km %.1f' % (lo, hi, a[thin & (dd >= lo) & (dd < hi)].sum()) for lo, hi in ((0, 150), (150, 300), (300, 500), (500, 1000), (1000, 9e9))) + ' Mkm2')
np.savez_compressed('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz', T=T.astype(np.float32), ty=ty, z=z, de=de.astype(np.float32), db=db.astype(np.float32))
