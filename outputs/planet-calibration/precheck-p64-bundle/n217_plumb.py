import time, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(217)
B = {'crust_kappa_m2_yr': 0., 'trench_polarity': 'pair', 'trench_flux': 'consumed', 'arc_H_ocean_m': 1000., 'arc_accretion': True,
     'collision_spread': 'earth', 'eclogite_rule': 'foundering', 'dyn_rule': 'measured'}
sim = P.PlateSim(g, seed=1012, params=B)
t0 = time.time()
for s in range(12):
    sim.step()
print('12 steps %.0f s; z max %.0f m; crust max %.1f km; spread %.2e m3; foundered %.2e m3' % (time.time()-t0, sim.elevation().max(),
      sim.world['thick'][sim.world['cont']].max()/1e3, sim.erosion_totals.get('collision_spread_m3', 0), sim.erosion_totals.get('eclogite_foundered_m3', 0)))
