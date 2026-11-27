import sys, json, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid.fields import nearest_nodes
from planet_earth.reference import ice_mask
R = 6371.0
g = CubeSphere(217); xyz = g.xyz; A = g.area
def dist_to_ocean(land):
    nb = g.nbr; okn = nb >= 0
    coast = (~land) & (okn & land[np.where(okn, nb, 0)]).any(1)
    oc = xyz[coast]; d = np.full(len(xyz), 0.0); idx = np.where(land)[0]
    for s in range(0, len(idx), 2000):
        c = xyz[idx[s:s+2000]] @ oc.T
        d[idx[s:s+2000]] = np.arccos(np.clip(c.max(1), -1, 1))*R
    return d
def wmed(v, w, q=50):
    o = np.argsort(v); cw = np.cumsum(w[o]); return float(v[o][np.searchsorted(cw, cw[-1]*q/100)])
BINS = [0, 50, 150, 300, 600, 1200, 1e9]
def profile(z, land, T=None, Tmask=None, keep=None):
    keep = np.ones(len(z), bool) if keep is None else keep
    d = dist_to_ocean(land); out = []
    tot = A[land & keep].sum()
    for a, b in zip(BINS[:-1], BINS[1:]):
        m = land & keep & (d >= a) & (d < b)
        if not m.any(): out.append(None); continue
        r = {'bin_km': [a, b], 'area_share': round(float(A[m].sum()/tot), 3), 'z_p50': round(wmed(z[m], A[m])),
             'z_p25': round(wmed(z[m], A[m], 25)), 'z_p75': round(wmed(z[m], A[m], 75)),
             'lt200': round(float(A[m & (z < 200)].sum()/A[m].sum()), 3),
             'gt1000': round(float(A[m & (z > 1000)].sum()/A[m].sum()), 3)}
        if T is not None:
            mt = m & Tmask
            r['T_p50_km'] = round(wmed(T[mt], A[mt])/1e3, 1) if mt.any() else None
        out.append(r)
    return out
lat = np.degrees(np.arcsin(xyz[:, 2])); lon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
# Earth
ze = np.load('outputs/planet-reference/earth-n217.npz')['elevation_m'].astype(float)
Tc = np.loadtxt('data/crust1/crsthk.xyz')[:, 2].reshape(180, 360)*1e3
i = np.clip((90-lat).astype(int), 0, 179); j = np.clip((lon+180).astype(int), 0, 359)
Te = Tc[i, j]
res = {'earth': profile(ze, ze > 0, Te, Te >= 20e3, ~ice_mask(g))}
# planet, plate layer at 0 Myr (46 km)
s = np.load('outputs/planet-plates/freeze-20261031/seed1001/snapshots/t+0000.0.npz')
zp = s['elevation_m'].astype(float)   # already relative to sea level
Tp = (s['crust_thickness_m'] + s['sediment_m']).astype(float)
res['planet_plates'] = profile(zp, zp > 0, Tp, Tp >= 20e3)
# planet landscape at 11.53 km, onto n = 217
L = np.load('outputs/planet-landscape/freeze-20261031-seed1001/present.npz')
g8 = CubeSphere(868); nn = nearest_nodes(g, g8.xyz)
zs = np.bincount(nn, L['elevation_m']*g8.area, len(xyz)) / np.bincount(nn, g8.area, len(xyz))
ocs = np.bincount(nn, L['ocean']*g8.area, len(xyz)) / np.bincount(nn, g8.area, len(xyz))
res['planet_landscape'] = profile(zs, ocs < 0.5, Tp, Tp >= 20e3)
for k, v in res.items():
    print(k)
    for r in v:
        print('  ', r)
json.dump(res, open('outputs/planet-calibration/precheck-p53-freeze/coast_profile.json', 'w'), indent=1)
