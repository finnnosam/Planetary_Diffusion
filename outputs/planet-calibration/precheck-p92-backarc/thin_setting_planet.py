"""Pass 92: the planet's inland thin crust by its setting, as thin_setting.py reads Earth's: thin = total crust +
sediment 20-30 km on the label-free continent (>= 20 km); of the thin crust more than 200 km from the continent's edge,
the share within 600 / 1,000 km of a subduction boundary node (either side), against all continent in that band.
Earth (Bird 2003 steps): 0.20 / 0.57 against 0.07 / 0.19. Six-seed means over the snapshots -40 ... 0 Myr. Args: arm
folders. Run from the repo root."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
g = CubeSphere(217)
SUB = 2
a = g.area
for arm in sys.argv[1:]:
    rows = []
    for sd in sorted(Path(arm).glob('seed*')):
        for f in sorted((sd/'snapshots').glob('t*.npz')):
            if f.name.endswith('-plates.json'):
                continue
            q = np.load(f)
            t = float(q['time_myr'])
            if t < -40.5:
                continue
            tot = (q['crust_thickness_m']+q['sediment_m'])/1e3
            whole = q['continental'].astype(bool) & (tot >= 20)
            d = gl.distance_to(~whole, g.nbr, g.dist)/1e3
            ds = gl.distance_to(q['boundary_type'] == SUB, g.nbr, g.dist)/1e3
            inl = whole & (d >= 200); thin = inl & (tot < 30)
            T, A = a[thin].sum(), a[inl].sum()
            rows.append([a[whole & (tot < 30)].sum()/a[whole].sum(), T/max(a[whole & (tot < 30)].sum(), 1),
                         a[thin & (ds < 600)].sum()/max(T, 1), a[thin & (ds < 1000)].sum()/max(T, 1),
                         a[inl & (ds < 600)].sum()/max(A, 1), a[inl & (ds < 1000)].sum()/max(A, 1)])
    if rows:
        m = np.mean(rows, 0)
        print(f'{Path(arm).name}: thin share {m[0]:.3f}; of it beyond 200 km {m[1]:.2f}; inland thin within 600 / 1,000 km '
              f'of subduction {m[2]:.2f} / {m[3]:.2f}; all inland continent {m[4]:.2f} / {m[5]:.2f}  ({len(rows)} snapshots)')
print('Earth: thin share 0.155; of it beyond 200 km 0.27; inland thin within 600 / 1,000 km 0.20 / 0.57; all inland 0.07 / 0.19')
