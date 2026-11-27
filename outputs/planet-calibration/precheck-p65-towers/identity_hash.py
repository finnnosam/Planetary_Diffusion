import sys, hashlib, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
g = CubeSphere(48)
out = []
for seed in (3, 11):
    sim = PlateSim(g, seed=seed, params={'threads': 1})
    for s in range(70):
        sim.step()
        st = sim.state()
        h = hashlib.md5()
        for k in sorted(st):
            v = np.asarray(st[k])
            if v.dtype.kind in 'fiub':
                h.update(np.ascontiguousarray(v).tobytes())
        out.append(h.hexdigest())
open(sys.argv[1], 'w').write('\n'.join(out))
print(len(out), 'steps')
