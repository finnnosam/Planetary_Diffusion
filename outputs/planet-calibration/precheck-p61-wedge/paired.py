# Pass 61 paired probe: the same state (seed 1010, combo-02's settings, 60 steps), then 20 steps with and without S.
import copy, json, sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_landscape import graph_landscape as gl
from planet_tectonics.isostasy import sediment_column
g = CubeSphere(217); a = g.area
base = {'crust_kappa_m2_yr': 0., 'trench_polarity': 'pair', 'trench_flux': 'consumed', 'arc_H_ocean_m': 1000.,
        'arc_accretion': True, 'lithosphere': 'none', 'dyn_rule': 'measured', 'threads': 1}
sim = P.PlateSim(g, seed=1010, params=base)
for _ in range(60):
    sim.step()
def read(s, t0):
    w = s.world; e = s.erosion_totals
    col = w['thick']+sediment_column(w['sed'], w['cont'], s.p)[0]
    cont = col >= 20e3; oc = ~cont
    sb = sediment_column(w['sed'], w['cont'], s.p)[0]
    dc = gl.distance_to(cont, g.nbr, g.dist, 1.1e8)
    bt = s.btype; dt = gl.distance_to(bt == 2, g.nbr, g.dist, 3e6)
    V = lambda m: float((sb[m]*a[m]).sum()/1e15)
    d = {k: (e.get(k, 0.)-t0.get(k, 0.))/1e9/20e6 for k in ('ocean_sediment_m3', 'deep_sea_sediment_m3', 'margin_wedge_m3',
                                                            'margin_wedge_spill_m3', 'sediment_subducted_m3', 'sediment_accreted_m3')}
    d.update({'ocean_sed_Mkm3': V(oc), 'far_share': V(oc & (dc >= 1e6))/V(oc), 'trench150_mean_km': V(oc & (dt < 150e3))/(a[oc & (dt < 150e3)].sum()/1e12),
              'by_dist': [V(oc & (dc >= lo*1e3) & (dc < hi*1e3)) for lo, hi in ((0, 200), (200, 500), (500, 1000), (1000, 1e6))]})
    return d
t0 = copy.deepcopy(sim.erosion_totals)
arm = sys.argv[1]
sim.p['margin_wedge'] = arm == 'S'                # the same state either way (deterministic on one thread)
for _ in range(20):
    sim.step()
r = read(sim, t0)
print(arm, json.dumps(r, default=float))
json.dump(r, open(f'/home/claude/pre/diag/paired_{arm}.json', 'w'), indent=1, default=float)
