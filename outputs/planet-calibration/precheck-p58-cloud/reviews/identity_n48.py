import sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
out, steps = sys.argv[1], int(sys.argv[2])
extra = eval(sys.argv[3]) if len(sys.argv) > 3 else {}
g = CubeSphere(48)
s = P.PlateSim(g, seed=1001, params={'threads': 1, **extra})
for i in range(steps):
    s.step()
st = s.state()
np.savez(out, **{k: np.asarray(v) for k, v in st.items() if isinstance(v, (np.ndarray, float, int, np.floating, np.integer))})
print('done', s.t, len(s.plates))
