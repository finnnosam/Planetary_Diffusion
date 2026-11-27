import sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(48)
sim = P.PlateSim(g, seed=1001, params={'threads': 1})
out = {}
for s in range(70):
    sim.step()
    st = sim.state()
    for k, v in st.items():
        if isinstance(v, np.ndarray):
            out[f'{s:02d}_{k}'] = v
np.savez(sys.argv[1], **out)
