"""Pass 54 diagnostic: control vs collapse flow off (crust_kappa_m2_yr 0), cloud seed 1007. Read-only.
Present day and the mean over the saved snapshots of the last 40 Myr (every 5 Myr)."""
import sys, json, glob, re
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import total_thickness_m, wpct
from planet_reviews.interior_probe import slope_m_km, edge_depth_km

g = CubeSphere(217); A = g.area; R = 6371.0
BINS = [0, 50, 150, 300, 600, 1200, 1e9]


def dist_to_ocean(land):
    nb = g.nbr; ok = nb >= 0
    coast = (~land) & (ok & land[np.where(ok, nb, 0)]).any(1)
    oc = g.xyz[coast]; d = np.zeros(g.size); idx = np.flatnonzero(land)
    for i in range(0, len(idx), 2000):
        d[idx[i:i+2000]] = np.arccos(np.clip((g.xyz[idx[i:i+2000]] @ oc.T).max(1), -1, 1))*R
    return d


def metrics(f, p):
    z = np.load(f); e = z['elevation_m'].astype(float); land = e > 0
    T = total_thickness_m(z, p); c = T >= 20e3
    d = dist_to_ocean(land)
    prof = [wpct(e[land & (d >= a) & (d < b)], A[land & (d >= a) & (d < b)], 50) for a, b in zip(BINS[:-1], BINS[1:])]
    far = land & (d >= 300)
    s = slope_m_km(g, e); hi = land & (e > 2000)
    thick = c & (T > 45e3)
    dep, pieces = edge_depth_km(g, thick)
    return {'cont': A[c].sum()/A.sum(), 'gt45': A[thick].sum()/A[c].sum(), 'gt55': A[c & (T > 55e3)].sum()/A[c].sum(),
            'T_p50_km': wpct(T[c], A[c], 50)/1e3, 'land_frac': A[land].sum()/A.sum(),
            'land_p25': wpct(e[land], A[land], 25), 'land_p50': wpct(e[land], A[land], 50), 'land_p90': wpct(e[land], A[land], 90),
            'hi2km_of_land': A[hi].sum()/A[land].sum(), 'slope_hi_p50': wpct(s[hi], A[hi], 50) if hi.any() else np.nan,
            'depth45_p50': wpct(dep[thick], A[thick], 50) if thick.any() else np.nan,
            'rise_0_600': prof[4]-prof[0], 'lt200_beyond300': A[far & (e < 200)].sum()/max(A[far].sum(), 1.),
            'sea_datum': float(z['sea_level_m']), 'profile': prof}


for run in sys.argv[1:]:
    sd = glob.glob(f'{run}/seed*')[0]
    st = json.load(open(f'{sd}/stats.json')); p = st['params_effective']
    fs = sorted(glob.glob(f'{sd}/snapshots/t*.npz'), key=lambda f: float(re.search(r't([+-][\d.]+)\.npz', f).group(1)))
    late = [f for f in fs if float(re.search(r't([+-][\d.]+)\.npz', f).group(1)) >= -40]
    ms = [metrics(f, p) for f in late]
    now = ms[-1]
    print(f'== {run}  crust volume change {st["continental_crust_volume_change"]:+.3f}, with sediment {st["continental_volume_change"]:+.3f}; '
          f'area ledger ' + ', '.join(f'{k} {v:+.3f}' for k, v in st['present']['area_ledger_frac'].items()))
    for k in now:
        if k == 'profile':
            print('   profile 0 Myr (land p50 by distance inland, m):', [round(x) for x in now[k]])
            continue
        print(f'   {k:15s} 0 Myr {now[k]:9.3f}   last-40 mean {np.nanmean([m[k] for m in ms]):9.3f}')
