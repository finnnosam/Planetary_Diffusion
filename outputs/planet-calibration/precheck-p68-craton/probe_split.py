"""Pass 68 probe: why a whole-continent plate stops splitting under split_rule 'earth' (seed 1001, init_site ocean).
Logs per step, for each plate with continent: Ac (Mkm^2), L of the candidate line (km), the rate prod*Ac/max(L,200),
and the successful rifts. Output only; the run is the model's own."""
import json, sys, time
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
from planet_tectonics import events as E
from planet_tectonics.breakup import production_per_mkm2

seed, site, steps, out = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]
g = CubeSphere(217)
sim = PlateSim(g, seed=seed, t0=-250., params={'init_site': site, 'threads': 1})
log = []
orig = type(sim)._split_margin_km


def wrapped(self, own, lab):
    L = orig(self, own, lab)
    A = self.g.area
    Ac = A[own & self.world['cont']].sum()/1e12
    prod = production_per_mkm2(self.t-self.t_start)
    n1 = int((own & (lab == 1)).sum()); n0 = int((own & (lab == 0)).sum())
    pid = int(np.bincount(self.owner[own]).argmax())
    log.append({'t': float(self.t), 'plate': pid, 'area': float(A[own].sum()/1e12), 'Ac': float(Ac), 'L': L,
                'prod': prod, 'rate': prod*Ac/max(L, E.SPLIT_L_MIN_KM) if L > 0 else 0., 'n0': n0, 'n1': n1})
    return L


type(sim)._split_margin_km = wrapped
sim.state()
t0 = time.time()
for k in range(steps):
    sim.step(); sim.state()
    if k % 10 == 0:
        print(seed, site, sim.t, len(sim.plates), round(time.time()-t0), flush=True)
json.dump({'log': log, 'events': [{kk: (float(v) if isinstance(v, (float, np.floating)) else v) for kk, v in e.items()
                                   if kk in ('t', 'kind', 'plate', 'new_plate', 'area_mkm2')} for e in sim.events]},
          open(out, 'w'))
