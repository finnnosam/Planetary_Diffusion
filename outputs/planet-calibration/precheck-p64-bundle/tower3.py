# Pass 64: the no-flow bundle at n = 96 (the tower test's setting), with and without R, 15 and 40 steps.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(96)
B = {'crust_kappa_m2_yr': 0., 'trench_polarity': 'pair', 'trench_flux': 'consumed', 'arc_H_ocean_m': 1000., 'arc_accretion': True}
for name, prm in (('control', {}), ('bundle', B), ('bundle + R', {**B, 'eclogite_rule': 'foundering'})):
    sim = P.PlateSim(g, seed=1001, params=prm)
    for s in range(40):
        sim.step()
        if s in (14, 39):
            w = sim.world; z = sim.elevation(); c = w['cont']
            print(f"{name:12s} step {s+1}: z max {z.max():6.0f} m, crust max {w['thick'][c].max()/1e3:5.1f} km, > 65 km {(c & (w['thick'] > 65e3)).sum()} nodes, "
                  f"sed max {w['sed'].max()/1e3:5.1f} km, foundered {sim.erosion_totals.get('eclogite_foundered_m3', 0)/1e9/(s+1)/1e6:.2f} km3/yr", flush=True)
