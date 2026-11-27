import sys, numpy as np, ducc0
E = ducc0.sht.experimental
cnt = {}
def wrap(name):
    f = getattr(E, name)
    def g(*a, **k):
        r1 = f(*a, **k); bad = 0
        for _ in range(4):
            bad += not np.array_equal(r1, f(*a, **k))
        key = (name, k.get('lmax'), k.get('nthreads')); c = cnt.setdefault(key, [0, 0]); c[0] += 1; c[1] += bad > 0
        return r1
    setattr(E, name, g)
wrap('adjoint_synthesis_general'); wrap('synthesis_general')
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
sim = PlateSim(CubeSphere(48), seed=1001, t0=-250., params={'threads': int(sys.argv[1])})
sim.state(); sim.run(10, lambda s: s.state(), 1)
for k, v in sorted(cnt.items(), key=str): print(k, f'{v[1]} of {v[0]} not repeatable')
