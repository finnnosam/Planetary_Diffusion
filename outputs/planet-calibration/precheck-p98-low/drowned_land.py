"""Pass 98 probe: the drowned continent's two parts (drowned_split.py) by distance from land (graph distance to the
nearest node above the sea): share of each part at 0-50 / 50-150 / 150-300 / > 300 km, and the part's room under the
graded shelf's ceiling (-(10 + 1.2 m/km d) within 150 km; none beyond) and under Earth's submerged-continent p50 by
distance from land (shelf_profile.py, ETOPO): mean m. Six seeds, 0 Myr. Read, not a rule. Usage: <arm> ..."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
g = CubeSphere(217); a = g.area
EB = np.array([0, 50, 100, 150, 200, 300, 400, 600, 1000.]); EP = np.array([-51, -147, -212, -236, -293, -692, -892, -1306.])
B = (0, 50, 150, 300, 1e9)
for arm in sys.argv[1:]:
    R = []
    for sd in sorted(Path(arm).glob('seed*')):
        q = np.load(sd/'snapshots'/'t+0000.0.npz')
        Tc = q['crust_thickness_m'].astype(float); S = q['sediment_bulk_m'].astype(float); z = q['elevation_m'].astype(float)
        c = (Tc+S) >= 20e3; dl = gl.distance_to(z >= 0, g.nbr, g.dist)/1e3
        ceil = np.where(dl <= 150, -(10+1.2*dl), -1e9)
        ear = EP[np.clip(np.searchsorted(EB, dl, 'right')-1, 0, len(EP)-1)]
        dr = c & (z >= -1000) & (z < -200); row = []
        for m in (dr & (Tc < 30e3), dr & (Tc >= 30e3)):
            W = a[m].sum()
            row += [a[m & (dl >= lo) & (dl < hi)].sum()/W for lo, hi in zip(B[:-1], B[1:])]
            row += [float((a[m]*np.maximum(ceil[m]-z[m], 0)).sum()/W), float((a[m]*np.maximum(ear[m]-z[m], 0)).sum()/W),
                    float(np.median(dl[m]))]
        R.append(row)
    M = np.mean(R, 0)
    print(Path(arm).name, '(6 seeds, 0 Myr): share at 0-50 / 50-150 / 150-300 / > 300 km from land; room under the graded ceiling / under Earth p50 (m); median km')
    for k, lab in ((0, 'drowned Tc < 30 '), (7, 'drowned Tc >= 30')):
        print(f'  {lab}: ' + ' / '.join(f'{v:.2f}' for v in M[k:k+4]) + f' | room {M[k+4]:.0f} / {M[k+5]:.0f} m | median {M[k+6]:.0f} km')
