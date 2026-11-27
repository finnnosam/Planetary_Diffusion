"""Run-to-run repeat check: n = 48, seed, threads T, K steps; md5 of every state() array per step (calls state() as simulate_plates does)."""
import sys, hashlib, json, numpy as np
n, seed, T, K, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
g = CubeSphere(n)
sim = PlateSim(g, seed=seed, t0=-250., params={'threads': T})
rows = []
def rec(s):
    st = s.state()
    rows.append({k: hashlib.md5(np.ascontiguousarray(v).tobytes()).hexdigest() for k, v in st.items() if isinstance(v, np.ndarray)})
    rows[-1]['_flowH'] = hashlib.md5(s.mantle.flow.H.tobytes()).hexdigest() if getattr(s, 'mantle', None) is not None and s.mantle.flow is not None else ''
rec(sim)
sim.run(K, rec, 1)
json.dump(rows, open(out, 'w'))
print('done', len(rows))
