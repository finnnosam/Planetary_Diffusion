import numpy as np, json
from planet_grid.sphere_graph import CubeSphere
g = CubeSphere(217); xyz = g.xyz; A = g.area; R = 6371.0
s = np.load('outputs/planet-plates/freeze-20261031/seed1001/snapshots/t+0000.0.npz')
sl = float(s['sea_level_m']); z = s['elevation_m'].astype(float); dyn = s['dynamic_topography_m'].astype(float)
Tc = s['crust_thickness_m'].astype(float); sed = s['sediment_m'].astype(float); sb = s['sediment_bulk_m'].astype(float)
land = z > 0   # elevation_m is relative to sea level; datum = z + sl
nb = g.nbr; okn = nb >= 0
coast = (~land) & (okn & land[np.where(okn, nb, 0)]).any(1)
oc = xyz[coast]; d = np.zeros(len(z)); idx = np.where(land)[0]
for i in range(0, len(idx), 2000):
    d[idx[i:i+2000]] = np.arccos(np.clip((xyz[idx[i:i+2000]] @ oc.T).max(1), -1, 1))*R
# Airy of crust alone, model's reference: 38 km -> 300 m (above the datum), rho 2830 / 3300
airy = 300 + (Tc - 38000)*(1 - 2830/3300)
print('sea level vs datum', round(sl), 'm')
def wm(v, w): o = np.argsort(v); c = np.cumsum(w[o]); return float(v[o][np.searchsorted(c, c[-1]/2)])
print('bin_km  z(sea)  z(datum)  airy(Tc,datum)  dyn  rest=z(datum)-airy-dyn  Tc_km  sed_km')
for a, b in zip([0,50,150,300,600,1200],[50,150,300,600,1200,1e9]):
    m = land & (d >= a) & (d < b); w = A[m]
    print(f'{a:>5}-{b:<6.0f}', *[round(wm(v[m], w)) for v in (z, z+sl, airy, dyn, z+sl-airy-dyn)], round(wm(Tc[m],w)/1e3,1), round(wm(sed[m],w)/1e3,2))
# least squares over land: z = c0 + c1*Tc + c2*dyn
M = land; X = np.stack([np.ones(M.sum()), Tc[M]/1e3, dyn[M], sed[M]/1e3], 1); wts = np.sqrt(A[M])
c = np.linalg.lstsq(X*wts[:,None], (z[M]+sl)*wts, rcond=None)[0]
print('fit z(datum) = %.0f + %.1f m/km Tc + %.2f dyn + %.0f m/km sed' % tuple(c))
