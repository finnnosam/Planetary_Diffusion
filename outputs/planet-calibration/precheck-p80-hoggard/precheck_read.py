"""Pass 80 (journal 104): read the cloud pre-check (seed 1011, 250 Myr; control and lith_frame 'lamb'): land
p25 / p50 / p90 at 0 Myr, land fraction, sea level, L by crust bin in the run's own frame, dL/dT, L's maximum over the
saved snapshots (bounded?), the highest node at 0 Myr. Read only, not pooled. Run from the repo root: precheck_read.py <run dirs>"""
import sys, glob, json, numpy as np
from planet_tectonics.lithosphere import thickness_km
from planet_reviews.crust_review import run_params
for run in sys.argv[1:]:
    sd = (glob.glob(f'{run}/seed1011')+glob.glob(f'{run}/*/seed1011'))[0]
    p = run_params(sd)
    q = np.load(f'{sd}/snapshots/t+0000.0.npz')
    z = q['elevation_m']; c = q['continental'].astype(bool); la = q['lithosphere_age_myr'].astype(float)
    T = (q['crust_thickness_m']+q['sediment_m'])/1e3
    land = z > 0
    m = c & np.isfinite(la) & (T >= 20); L = thickness_km(la[m], p)
    Lmax = max(np.nanmax(thickness_km(np.load(f)['lithosphere_age_myr'].astype(float), p))
               for f in sorted(glob.glob(f'{sd}/snapshots/t*.npz')))
    print(f'{run} (frame {p.get("lith_frame", "artemieva")}): land {land.mean():.3f}, p25/p50/p90 '
          f'{np.round(np.percentile(z[land], [25, 50, 90]))} m, highest node {z.max():.0f} m, L max over the saves {Lmax:.1f} km')
    print(f'  dL/dT {np.polyfit(T[m], L, 1)[0]:.2f} km/km, corr {np.corrcoef(T[m], L)[0, 1]:.2f}; L mean by crust bin: ' +
          ', '.join(f'{a}-{b} {L[(T[m] >= a) & (T[m] < b)].mean():.0f}' for a, b in
                    [(20, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)] if ((T[m] >= a) & (T[m] < b)).sum() > 20))
    print('  z p50 by crust bin: ' + ', '.join(f'{a}-{b} {np.median(z[c & (T >= a) & (T < b)]):.0f}' for a, b in
          [(20, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)] if (c & (T >= a) & (T < b)).sum() > 20))
