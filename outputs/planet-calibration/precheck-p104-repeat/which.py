import sys, numpy as np
import planet_tectonics.mantle_flow as mf
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
cnt = {}
def wrap(name):
    f = getattr(mf, name)
    def g(*a, **k):
        r1 = f(*a, **k); bad = 0
        for _ in range(8):
            r2 = f(*a, **k); bad += not np.array_equal(r1, r2)
        c = cnt.setdefault(name, [0, 0]); c[0] += 1; c[1] += bad > 0
        return r1
    setattr(mf, name, g)
for nm in ('point_alm', 'synthesis', 'gradient'): wrap(nm)
sim = PlateSim(CubeSphere(48), seed=1001, t0=-250., params={'threads': int(sys.argv[1])})
sim.state(); sim.run(12, lambda s: s.state(), 1)
print('threads', sys.argv[1], {k: f'{v[1]} of {v[0]} calls not repeatable' for k, v in cnt.items()})
