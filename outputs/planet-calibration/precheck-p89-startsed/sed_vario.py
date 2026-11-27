"""Correlation of normal scores vs separation: Earth's stable-crust sediment against random_field(lmax, slope)."""
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_landscape import graph_landscape as gl
from planet_reviews.interior_probe import GROUPS
from planet_earth.bird import earth_boundaries
from planet_earth.crust1 import C1, crust1_cell
from planet_tectonics.lithosphere import _nscore
g = CubeSphere(217); a = g.area; cell = crust1_cell(g)
T = np.loadtxt(C1/'crsthk.xyz')[:, 2][cell]*1e3; S = np.loadtxt(C1/'sedthk.xyz')[:, 2][cell]*1e3
ty = np.array([l.split() for l in (C1/'CNtype1-1.txt').read_text().splitlines() if l.strip()]).ravel()[cell]
dbnd = gl.distance_to(earth_boundaries(g), g.nbr, g.dist, 5000e3)/1e3
m = (T >= 20e3) & np.isin(ty, GROUPS['Precambrian/platform']) & (dbnd > 500)
idx = np.flatnonzero(m); rng = np.random.default_rng(0)
sub = rng.choice(idx, 4000, replace=False)
X = g.xyz[sub]; D = np.degrees(np.arccos(np.clip(X@X.T, -1, 1)))*111.2
bins = np.array([0, 100, 200, 400, 600, 1000, 1500, 2500])
def vario(z):
    zz = z[:, None]*z[None, :]
    return [float(zz[(D > lo) & (D <= hi)].mean()) for lo, hi in zip(bins[:-1], bins[1:])]
zs = _nscore(S[m], a[m]); zmap = np.zeros(g.size); zmap[m] = zs
print('Earth sed  ', np.round(vario(zmap[sub]), 2))
zc = np.zeros(g.size); zc[m] = _nscore((T-S)[m], a[m]); print('Earth crust', np.round(vario(zc[sub]), 2))
for lmax, sl in ((24, 2.), (12, 2.5), (40, 2.), (24, 1.)):
    r = []
    for s in range(3):
        f = fields.random_field(g.xyz, lmax, sl, np.random.default_rng(s)); r.append(vario(f[sub]))
    print(f'field {lmax},{sl}', np.round(np.mean(r, 0), 2))
print('---')
E = np.array(vario(zmap[sub]))
for sl in (.5, 1., 1.5, 2.):
    r = np.mean([vario(fields.random_field(g.xyz, 24, sl, np.random.default_rng(10+s))[sub]) for s in range(6)], 0)
    print(f'field 24,{sl}', np.round(r, 2), 'sse', round(float(((r-E)**2).sum()), 3))
print('--- normalised')
def corr(z):
    z = (z-z.mean())/z.std(); return np.array(vario(z))
E = corr(zmap[sub]); print('Earth', np.round(E, 2))
for lmax in (24, 60, 100):
    for sl in (.5, 1., 1.5):
        r = np.mean([corr(fields.random_field(g.xyz, lmax, sl, np.random.default_rng(10+s))[sub]) for s in range(6)], 0)
        print(f'field {lmax},{sl}', np.round(r, 2), 'sse', round(float(((r-E)**2).sum()), 3))
