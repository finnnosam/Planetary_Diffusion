import hashlib, sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
import planet_tectonics.sim as P
g = CubeSphere(48)
sim = P.PlateSim(g, seed=1001, params={'threads': 1})
h = hashlib.md5()
for i in range(70):
    sim.step()
    st = sim.state()
    for k in sorted(st):
        v = st[k]
        if isinstance(v, np.ndarray): h.update(k.encode()); h.update(np.ascontiguousarray(v).tobytes())
    if i % 10 == 9: print(i+1, h.hexdigest(), flush=True)
