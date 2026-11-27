"""Pass 101: the planet's counterpart of earth_pb.py: corr(field by degree band, distance to subduction nodes
(boundary_type SUBDUCTION)) at 0 Myr, on the whole air-loaded field (dynamic_topography_air_m, swells included) by
quadrature and band (so flow_placed 'earth_corr' reads its placed and plate-blind parts together), six seeds.
Read. Usage (repo root): python3 planet_pb.py <arm> ..."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_tectonics import plates as PL
from degree_split import corr, dist_to
g = CubeSphere(217); a = g.area
for arm in map(Path, sys.argv[1:]):
    R = []
    for sd in range(1001, 1007):
        q = np.load(arm/f'seed{sd}'/'snapshots'/'t+0000.0.npz')
        air = q['dynamic_topography_air_m'].astype(float); lm = 31; L = mf.alm_l(lm)
        H = mf.point_alm(g.xyz, air*a/g.radius**2, lm)
        d = dist_to(g, np.flatnonzero(q['boundary_type'] == PL.SUBDUCTION))
        row = [corr(mf.synthesis(np.where((L >= lo) & (L <= hi), H, 0), lm, g.xyz), d, a) for lo, hi in ((1, 3), (4, 10), (11, lm))]
        row.append(corr(air, d, a))
        R.append(row)
    M = np.mean(R, 0)
    print(f'{arm.name} (6 seeds, 0 Myr): corr with distance to subduction, air field l 1-3 / 4-10 / 11-31 {M[0]:+.2f} / {M[1]:+.2f} / {M[2]:+.2f}; whole air field {M[3]:+.2f}')
