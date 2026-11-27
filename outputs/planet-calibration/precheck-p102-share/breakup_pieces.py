"""Pass 102 follow-up probe: how far the supercontinent broke up. Per arm and seed at 0 Myr: the largest plate's area
and its continental area (label-free, crust + bulk sediment >= 20 km), and the continent's connected pieces (on the
graph's neighbours) >= 1 Mkm^2: count and the three largest (Mkm^2). Earth (CRUST1.0 on the same grid, if present)
beside it. Read, not a rule. Usage (repo root): python3 breakup_pieces.py <arm> ..."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.graphops import components
g = CubeSphere(217); A = g.area/1e12
for arm in sys.argv[1:]:
    print(Path(arm).name)
    for s in range(1001, 1007):
        q = np.load(Path(arm)/f'seed{s}'/'snapshots'/'t+0000.0.npz')
        pid = q['plate_id']; c = (q['crust_thickness_m']+q['sediment_bulk_m']) >= 20e3
        ids, inv = np.unique(pid, return_inverse=True)
        pa = np.bincount(inv, weights=A); big = np.argmax(pa)
        lab = components(c, g.nbr)
        pa_c = np.bincount(lab[c], weights=A[c])
        pieces = np.sort(pa_c[pa_c >= 1])[::-1]
        print(f'  seed{s}: largest plate {pa[big]:5.0f} Mkm2, its continent {A[(inv == big) & c].sum():5.0f}; continent '
              f'{A[c].sum():4.0f} in {pieces.size} pieces >= 1 Mkm2, largest {", ".join(f"{x:.0f}" for x in pieces[:3])}')
