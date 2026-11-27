# Pass 64: which term thickens the tallest columns (n = 96, bundle + Earth spread): per step, the volume added to
# continental crust thicker than 70 km by subduction thickening (+ volcanic), by collision volume, and in total.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(96); a = g.area
B = {'crust_kappa_m2_yr': 0., 'trench_polarity': 'pair', 'trench_flux': 'consumed', 'arc_H_ocean_m': 1000., 'arc_accretion': True, 'collision_spread': 'earth'}
sim = P.PlateSim(g, seed=1001, params=B)
rec = {}
orig_sub = sim._subduction_thickening
def sub_wrap():
    out = orig_sub(); rec['sub'] = out.copy(); return out
sim._subduction_thickening = sub_wrap
orig_sp = sim._spread_collision
def sp_wrap(dV):
    out = orig_sp(dV); rec['coll'] = out/a; return out
sim._spread_collision = sp_wrap
for s in range(40):
    rec.clear()
    th0 = sim.world['thick'].copy()
    sim.step()
    w = sim.world; thick = w['cont'] & (w['thick'] > 70e3)
    S = rec.get('sub', np.zeros(g.size)); C = rec.get('coll', np.zeros(g.size))
    if s % 4 == 3:
        i = int(np.argmax(np.where(w['cont'], w['thick'], 0)))
        print(f"step {s+1}: > 70 km {thick.sum()} nodes; on them sub {S[thick].mean():6.0f} m, coll {C[thick].mean():6.0f} m; max node: T {w['thick'][i]/1e3:.0f} km, sub {S[i]:.0f}, coll {C[i]:.0f}, btype {sim.btype[i]}", flush=True)
