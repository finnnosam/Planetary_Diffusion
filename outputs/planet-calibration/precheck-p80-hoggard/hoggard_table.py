"""Pass 80 (journal 104): Earth's seismic lithosphere by tectonothermal age class, as a monotone table (C's 'lamb'
frame). Class means (area-weighted, n = 217 nodes on lithosphere.earth_continent's domain) of each Hoggard et al.
2020 model, pooled where an older class is thinner than a younger one (pool-adjacent-violators, area-weighted:
conductive cooling never thins a lithosphere as its thermal age grows; pooling keeps the area mean). Each on the whole
continent and on crust <= 45 km (unthickened: the planet thickens collided lithosphere itself, 'thicken'); the
reference is the table's area mean over the whole continent, as C's.
Run from the repo root after hoggard_law.py (imports it)."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(__file__))
import hoggard_law as h


def pava(x, w):
    blocks = [[x[i], w[i], [i]] for i in range(len(x))]
    i = 0
    while i < len(blocks)-1:
        if blocks[i][0] > blocks[i+1][0]+1e-12:
            a, b = blocks[i], blocks[i+1]
            blocks[i] = [(a[0]*a[1]+b[0]*b[1])/(a[1]+b[1]), a[1]+b[1], a[2]+b[2]]; del blocks[i+1]; i = max(i-1, 0)
        else:
            i += 1
    out = np.empty(len(x))
    for v, _, idx in blocks:
        out[idx] = v
    return out


for model, dom in [(mo, d) for mo in h.MODELS for d in ('all', 'T<=45')]:
    L = h.sample(model); m0 = h.cont & np.isfinite(L) & (L <= 350)
    m = m0 & (h.T <= 45) if dom == 'T<=45' else m0
    ages = np.array(sorted(set(h.age[m]))); mu = []; wa = []
    for t0 in ages:
        k = m & (h.age == t0); mu.append((L[k]*h.a[k]).sum()/h.a[k].sum()); wa.append(h.a[k].sum())
    mu = np.array(mu); wa = np.array(wa); iso = pava(mu, wa)
    k = m & (h.age == ages[0]) & (h.T <= 45)
    lawall = np.interp(h.age[m0], ages, iso); ref_all = (lawall*h.a[m0]).sum()/h.a[m0].sum()
    print(f'{model:12s} [{dom}] ages {ages.astype(int).tolist()}  ref on the whole continent {ref_all:.2f}')
    print(f'{"":12s} means {np.round(mu, 1).tolist()}')
    print(f'{"":12s} monotone {np.round(iso, 2).tolist()}  ref {(iso*wa).sum()/wa.sum():.2f}'
          f'  youngest class without crust > 45 km: {(L[k]*h.a[k]).sum()/h.a[k].sum():.1f}')
