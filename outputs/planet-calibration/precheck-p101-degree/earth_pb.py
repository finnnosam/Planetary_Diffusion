"""Pass 101: Earth's counterpart of the design check P-b: the area-weighted correlation of Earth's residual topography
(Holdt 2022 degree 1-10, 1-30; VD2020's posterior mean by degree band, the reader of earth_degree_corr.py) with the
distance to subduction (Bird 2003 SUB step midpoints), beside the placed share of Earth's power under flow_placed
'earth_corr'. Read. Usage (repo root): python3 earth_pb.py"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.trench_origin import read_dbf, _xyz
from planet_tectonics import mantle_flow as mf
from planet_tectonics.mantle import measured_power_air_m2, placed_share
from degree_split import holdt, corr
from earth_degree_corr import vd2020_alm, LMAX
g = CubeSphere(217); a = g.area
s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf'); sub = np.asarray(s['STEPCLASS']) == 'SUB'
m = _xyz(np.asarray(s['STARTLONG'])[sub], np.asarray(s['STARTLAT'])[sub])+_xyz(np.asarray(s['FINALLONG'])[sub], np.asarray(s['FINALLAT'])[sub])
m /= np.linalg.norm(m, axis=1, keepdims=True)
d = np.empty(g.size)
for i in range(0, g.size, 4000):
    d[i:i+4000] = np.arccos(np.clip((g.xyz[i:i+4000]@m.T).max(1), -1, 1))*g.radius/1e3
h10, h30 = holdt(g, 'degree_1_to_10.grd'), holdt(g, 'degree_1_to_30.grd')
print(f'Holdt 1-10 / 11-30 / 1-30: corr with distance to subduction {corr(h10, d, a):+.2f} / {corr(h30-h10, d, a):+.2f} / {corr(h30, d, a):+.2f}')
V = vd2020_alm(False, True); L = mf.alm_l(LMAX)
for lo, hi in ((1, 3), (4, 10), (11, 30), (1, 30)):
    f = mf.synthesis(np.where((L >= lo) & (L <= hi), V, 0), LMAX, g.xyz)
    print(f'VD2020 mean l {lo}-{hi}: corr with distance to subduction {corr(f, d, a):+.2f}')
pw = measured_power_air_m2(); sh = placed_share(30)
print(f'placed share of the measured power: all {(pw*sh).sum()/pw.sum():.3f}; l 1-3 {(pw[1:4]*sh[1:4]).sum()/pw[1:4].sum():.3f}; '
      f'l 4-10 {(pw[4:11]*sh[4:11]).sum()/pw[4:11].sum():.3f}; l 11-30 {(pw[11:]*sh[11:]).sum()/pw[11:].sum():.3f}')
