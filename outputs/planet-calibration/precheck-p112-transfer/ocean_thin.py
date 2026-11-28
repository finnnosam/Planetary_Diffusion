"""Pass 112 probe (from precheck-p111-confirm/ocean_thin.py, seeds read from the folder): oceanic crust thinner than
ocean_thick_m (7 km) in every saved snapshot, which no ledger books.
Usage (repo root): python ocean_thin.py <arm folder> [<arm folder> ...]
Per arm and seed: the thinnest oceanic node over all snapshots, the snapshots with any node below 6,999 m, and the
largest deficit below 7 km (km^3)."""
import sys
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere

g = CubeSphere(217)
a = g.area
for arm in sys.argv[1:]:
    print(Path(arm).name)
    for sd in sorted(Path(arm).glob('seed*')):
        mins, bad, worst = [], 0, 0.
        snaps = sorted((sd/'snapshots').glob('t*.npz'))
        for f in snaps:
            with np.load(f) as d:
                oc = ~d['continental']
                T = d['crust_thickness_m'][oc].astype(float)
            mins.append(T.min())
            if (T < 6999).any():
                bad += 1
                worst = max(worst, float((np.clip(7000-T, 0, None)*a[oc]).sum())/1e9)
        print(f'  {sd.name}: {len(snaps)} snapshots, thinnest oceanic node {min(mins):.0f} m, snapshots below 6,999 m '
              f'{bad}, largest deficit {worst:.1f} km3')
