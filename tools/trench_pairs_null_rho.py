"""Forty-fourth pass: share of boundary length with rho < 0.5 on Earth's plates under random pole directions (each plate
keeps its NNR rate); the null for planet_reviews.trench_pairs. Run from the repo root with PYTHONPATH set."""
import numpy as np
from planet_reviews import trench_pairs as tp, force_balance_earth as fbe
g, pid, names, _ = fbe.earth_plates(217)
w = fbe.earth_omega(g, pid, names)
st = fbe.earth_steps(names, fbe._rings())
E = tp.earth_segments()
x = st["x"]*fbe.R*1e3*fbe.YR; L = E["L"]
cls = E['cls']; oo = E['io'] | E['ocb']
mags = [np.linalg.norm(w[k]) for k in w]
rng = np.random.default_rng(1)
res = []
for r in range(40):
    wr = {k: (lambda v: v/np.linalg.norm(v)*np.linalg.norm(w[k]))(rng.normal(size=3)) for k in w}
    wr = fbe.remove_net_rotation(g, pid, wr)
    va = np.cross(np.stack([wr[int(k)] for k in st['A']]), x); vb = np.cross(np.stack([wr[int(k)] for k in st['B']]), x)
    rho = tp.rho(va, vb)
    row = [tp.share(rho[oo] < .5, L[oo])]
    for c in ('OSR', 'CRB', 'OTF', 'CTF', 'SUB', 'OCB', 'CCB'):
        m = cls == c; row.append(tp.share(rho[m] < .5, L[m]))
    res.append(row)
res = np.array(res)
print('Earth plates, random pole directions, each plate keeps its NNR rate (40 draws): share of length at rho < 0.5, mean (p5-p95)')
for name, col in zip(['oo convergence', 'OSR', 'CRB', 'OTF', 'CTF', 'SUB', 'OCB', 'CCB'], res.T):
    print(f'  {name:15s} {col.mean():.2f} ({np.percentile(col, 5):.2f}-{np.percentile(col, 95):.2f})')
