import numpy as np, sys
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(96)
OLD = {'crust_kappa_m2_yr': 1000., 'trench_polarity': 'node', 'trench_flux': 'band', 'arc_H_ocean_m': 300., 'arc_accretion': False}
arms = {'new': {}, 'old': OLD}
for k in OLD: arms['old but '+k+' new'] = {kk: v for kk, v in OLD.items() if kk != k}
for name, prm in arms.items():
    sim = P.PlateSim(g, seed=1001, params=prm)
    for _ in range(15): sim.step()
    w = sim.world; i = int(np.argmax(w['sed']))
    print(f"{name:40s} sed max {w['sed'].max()/1e3:6.1f} km at node {i}: cont {bool(w['cont'][i])}, crust {w['thick'][i]/1e3:.1f} km, sed>10km nodes {(w['sed']>10e3).sum()}", flush=True)
