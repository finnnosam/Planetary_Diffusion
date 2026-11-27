"""Pass 102 follow-up probe (the user saw seed 1001's supercontinent stay whole under 'earth_corr'): the most continent
(label-free, crust + bulk sediment >= 20 km) held by one plate, by seed, at -200 / -150 / -100 / -50 / 0 Myr where saved
(sweeps save the last 40 Myr: earlier times read the stats series' largest plate instead). Earth: Bird 2003's plates
on CRUST1.0 (42 / 35 / 27 Mkm^2 on the three largest). Read, not a rule. Usage: python3 plate_continent.py <arm> ..."""
import sys, json
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
A = CubeSphere(217).area/1e12
for arm in sys.argv[1:]:
    out, big = [], []
    for s in range(1001, 1007):
        d = np.load(Path(arm)/f'seed{s}'/'snapshots'/'t+0000.0.npz')
        c = (d['crust_thickness_m']+d['sediment_bulk_m']) >= 20e3
        ids, inv = np.unique(d['plate_id'], return_inverse=True)
        out.append(np.bincount(inv, weights=A*c).max())
        ser = json.load(open(Path(arm)/f'seed{s}'/'stats.json'))['series']
        big.append([round(max(ser[i]['plate_areas_mkm2'])) for i in (50, 100, 150, 200)])
    print(Path(arm).name, '| continent on one plate at 0 Myr by seed', [round(x) for x in out], 'mean', round(np.mean(out)),
          '| largest plate at -200/-150/-100/-50 Myr by seed', big)
