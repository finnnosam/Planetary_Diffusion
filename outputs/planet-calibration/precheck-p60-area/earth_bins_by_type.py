# Pass 60: Earth's crust (ice out) by thickness bin and CRUST1.0 type group, Mkm2; and z p50 per group in the bin
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_earth.reference import ice_mask
g = CubeSphere(217); a = g.area/1e12; ice = ice_mask(g)
e = np.load('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz'); T, ty, z, de = e['T'], e['ty'], e['z'], e['de']
G = {'shelf': ['C-'], 'slope/trans': ['S-', 'T-', 'U-'], 'extended': ['M-', 'N-', 'X-', 'V2', 'Y1', 'Y2'], 'arc': ['J-', 'K-', 'L1', 'L2'],
     'ocean plat': ['B-', 'W-', 'A0', 'A1', 'V1', 'Y3'], 'orogen': ['O-', 'P-', 'Q-', 'R1', 'R2'],
     'Phanero': ['Z1', 'Z2'], 'platform': ['D-', 'E-'], 'Proter': ['H1', 'H2', 'I1', 'I2'], 'Archean': ['F-', 'G1', 'G2']}
B = [(20, 25), (25, 27.5), (27.5, 30), (30, 32.5), (32.5, 35), (35, 40), (40, 45), (45, 99)]
print('%-12s' % 'bin km' + ''.join('%11s' % k for k in G) + '%9s' % 'all')
for lo, hi in B:
    m = ~ice & (T >= lo*1e3) & (T < hi*1e3)
    print('%5.1f-%-6.1f' % (lo, hi) + ''.join('%11.1f' % a[m & np.isin(ty, v)].sum() for v in G.values()) + '%9.1f' % a[m].sum())
m = ~ice & (T >= 27.5e3) & (T < 32.5e3)
print('27.5-32.5 km: share within 150 km of a continental edge %.2f, submerged %.2f' % (a[m & (de < 150)].sum()/a[m].sum(), a[m & (z < 0)].sum()/a[m].sum()))
