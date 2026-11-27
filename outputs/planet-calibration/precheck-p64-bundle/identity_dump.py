# Pass 64: default path step by step (n = 48, 70 steps, threads 1); writes a hash per step for comparing two trees.
import sys, hashlib, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(48); sim = P.PlateSim(g, seed=5, params={'threads': 1})
out = []
for i in range(70):
    sim.step()
    h = hashlib.md5()
    for k in sorted(sim.world): h.update(np.ascontiguousarray(sim.world[k]).tobytes())
    h.update(sim.owner.tobytes()); out.append(h.hexdigest())
open(sys.argv[1], 'w').write('\n'.join(out)+'\n')
