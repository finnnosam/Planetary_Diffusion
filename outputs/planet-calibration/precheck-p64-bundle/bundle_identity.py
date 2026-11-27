# Pass 64: a bare PlateSim (the new PARAMS) against the old PARAMS with the bundle passed explicitly, step by step
# (n = 48, 70 steps, threads 1): identical means nothing reads the old defaults behind PARAMS.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
B = {'crust_kappa_m2_yr': 0., 'trench_polarity': 'pair', 'trench_flux': 'consumed', 'arc_H_ocean_m': 1000., 'arc_accretion': True}
g = CubeSphere(48)
a = P.PlateSim(g, seed=5, params={'threads': 1})
b = P.PlateSim(g, seed=5, params={'threads': 1, **B})
assert all(a.p[k] == v for k, v in B.items())
for i in range(70):
    a.step(); b.step()
    for k in a.world:
        x, y = np.asarray(a.world[k]), np.asarray(b.world[k])
        if not np.array_equal(x, y, equal_nan=True):
            raise SystemExit(f'step {i}: world[{k}] differs')
    assert np.array_equal(a.owner, b.owner)
print('identical over 70 steps; plates', len(a.plates))
