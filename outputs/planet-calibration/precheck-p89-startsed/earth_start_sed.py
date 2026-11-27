"""Earth's sediment on the start's crust domain (CRUST1.0 Precambrian/platform, >500 km from Bird's boundaries,
crust >= 20 km), n = 217: bulk and solid quantiles, rank correlation with crystalline crust, by distance from the
label-free edge; and on the whole label-free continent for comparison."""
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import wpct
from planet_reviews.interior_probe import GROUPS
from planet_earth.bird import earth_boundaries
from planet_earth.crust1 import C1, crust1_cell
from planet_tectonics.isostasy import _solid_thickness
from planet_tectonics.lithosphere import _nscore
import planet_tectonics.plates as P
g = CubeSphere(217); a = g.area; cell = crust1_cell(g)
T = np.loadtxt(C1/'crsthk.xyz')[:, 2][cell]*1e3; S = np.loadtxt(C1/'sedthk.xyz')[:, 2][cell]*1e3
ty = np.array([l.split() for l in (C1/'CNtype1-1.txt').read_text().splitlines() if l.strip()]).ravel()[cell]
dbnd = gl.distance_to(earth_boundaries(g), g.nbr, g.dist, 5000e3)/1e3
m = (T >= 20e3) & np.isin(ty, GROUPS['Precambrian/platform']) & (dbnd > 500)
lf = T >= 20e3
dco = gl.distance_to(~lf, g.nbr, g.dist, 5000e3)/1e3
p = P.PARAMS
Ss = _solid_thickness(S, p['phi0_continent'], p['compaction_c_continent_per_km'])
Q = (.01,.05,.10,.15,.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.80,.85,.90,.95,.99)
print('stable area Mkm2', a[m].sum()/1e12, ' label-free', a[lf].sum()/1e12)
print('bulk quantiles (m):', [round(float(wpct(S[m], a[m], 100*q))) for q in Q])
print('solid quantiles (m):', [round(float(wpct(Ss[m], a[m], 100*q))) for q in Q])
print('mean bulk / solid m', (S[m]*a[m]).sum()/a[m].sum(), (Ss[m]*a[m]).sum()/a[m].sum())
print('label-free mean bulk', (S[lf]*a[lf]).sum()/a[lf].sum(), 'vol Mkm3', (S[lf]*a[lf]).sum()/1e15)
zs, zt = _nscore(S[m], a[m]), _nscore((T-S)[m], a[m])
print('rank corr (sed, crystalline) stable:', np.sum(a[m]*zs*zt)/np.sqrt(np.sum(a[m]*zs**2)*np.sum(a[m]*zt**2)))
for lo, hi in ((0,200),(200,500),(500,1000),(1000,1500),(1500,5000)):
    k = m & (dco >= lo) & (dco < hi)
    if k.sum() > 20:
        print(f'dist {lo}-{hi} km: share {a[k].sum()/a[m].sum():.3f} mean bulk {(S[k]*a[k]).sum()/a[k].sum():.0f} p50 {wpct(S[k],a[k],50):.0f}')
