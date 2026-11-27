"""Step-by-step identity: n = 48, seed 5, threads 1, 70 steps, state() before the first step and after every step.
Prints one hash per step of every state array. Usage: identity.py <repo> <out> [json params]."""
import sys, json, hashlib, importlib
repo, out = sys.argv[1], sys.argv[2]
extra = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
import glob
sys.path[:0] = [repo]+sorted(glob.glob(repo+'/packages/*/src'))
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(48)
sim = P.PlateSim(g, seed=5, params={'threads': 1, **extra})
def h(st):
    m = hashlib.md5()
    for k in sorted(st):
        v = st[k]
        if isinstance(v, np.ndarray):
            m.update(k.encode()); m.update(np.ascontiguousarray(v).tobytes())
        else:
            m.update(f'{k}={v!r}'.encode())
    return m.hexdigest()
lines = [h(sim.state())]
for i in range(70):
    sim.step(); lines.append(h(sim.state()))
open(out, 'w').write('\n'.join(lines)+'\n')
print(out, lines[-1])
