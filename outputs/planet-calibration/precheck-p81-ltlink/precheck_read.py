"""Pass 81 (journal 105): read the cloud pre-check (seed 1012, 250 Myr; control and start_lithosphere 'crust' +
belt_lithosphere 'shorten'): land fraction and p25 / p50 / p90 at 0 Myr, sea level, the highest node over the saves,
L's maximum over the saves (bounded?), L by crust bin (all / old >= 250 Myr / young), young lithosphere under crust
> 45 km, dL/dT, z p50 by crust bin, crust > 45 km share. Read only, not pooled. Run from the repo root:
precheck_read.py <run dirs>"""
import sys, glob, numpy as np
from planet_tectonics.lithosphere import thickness_km
from planet_reviews.crust_review import run_params
from planet_grid.sphere_graph import CubeSphere
A = CubeSphere(217).area
BINS = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]
for run in sys.argv[1:]:
    sd = (glob.glob(f'{run}/seed1012')+glob.glob(f'{run}/*/seed1012'))[0]
    p = run_params(sd)
    snaps = sorted(glob.glob(f'{sd}/snapshots/t*.npz'))
    zmax = Lmax = -1e9
    for f in snaps:
        s = np.load(f); zmax = max(zmax, float(s['elevation_m'].max()))
        Lmax = max(Lmax, float(np.nanmax(thickness_km(s['lithosphere_age_myr'].astype(float), p))))
    for snap in ('t-0250.0', 't+0000.0'):
        q = np.load(f'{sd}/snapshots/{snap}.npz')
        z = q['elevation_m']; c = q['continental'].astype(bool); la = q['lithosphere_age_myr'].astype(float)
        T = (q['crust_thickness_m']+q['sediment_m'])/1e3
        land = z > 0
        m = c & np.isfinite(la) & (T >= 20); L = thickness_km(la[m], p); Tm = T[m]; ag = la[m]; w = A[m]
        wm = lambda k: (L[k]*w[k]).sum()/max(w[k].sum(), 1e-30)
        print(f'{run} {snap} (start {p["start_lithosphere"]}, belt {p["belt_lithosphere"]}): land {land.mean():.3f}, '
              f'p25/p50/p90 {np.round(np.percentile(z[land], [25, 50, 90]))} m, sea {q["sea_level_m"]:.0f} m; crust > 45 km '
              f'{w[Tm > 45].sum()/w.sum():.3f}' + (f'; over the saves: highest node {zmax:.0f} m, L max {Lmax:.1f} km' if snap[1] == '+' else ''))
        print(f'  dL/dT {np.polyfit(Tm, L, 1)[0]:.2f} km/km, corr {np.corrcoef(Tm, L)[0, 1]:.2f}; L mean all / old / young by crust bin: ' +
              '; '.join(f'{a}-{b} {wm((Tm >= a) & (Tm < b)):.0f}/{wm((Tm >= a) & (Tm < b) & (ag >= 250)):.0f}/{wm((Tm >= a) & (Tm < b) & (ag < 250)):.0f}'
                        for a, b in BINS if ((Tm >= a) & (Tm < b)).sum() > 20))
        print('  z p50 by crust bin: ' + ', '.join(f'{a}-{b} {np.median(z[c & (T >= a) & (T < b)]):.0f}' for a, b in BINS
                                                     if (c & (T >= a) & (T < b)).sum() > 20))
        yt = (Tm > 45) & (ag < 250)
        print(f'  young lithosphere under crust > 45 km: area share of thick {w[yt].sum()/max(w[Tm > 45].sum(), 1e-30):.2f}, L {wm(yt):.0f} km')
