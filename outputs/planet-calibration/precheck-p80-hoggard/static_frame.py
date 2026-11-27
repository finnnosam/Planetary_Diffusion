"""Pass 80 (journal 104): static estimate of C's 'lamb' frame on a swept arm's 0 Myr states (six seeds): the change of
C's term on every continental node (lamb minus the run's own frame, same thermal ages; water-loaded below sea level),
by crust bin, and land p25 / p50 / p90 with the change added (no sea-level response, no trajectory: an estimate, not
a reading). Also the term's area mean over the continent in each frame (a non-zero mean moves the sea).
Run from the repo root: static_frame.py <arm folder>"""
import sys, numpy as np
from planet_tectonics.lithosphere import elevation_term, thickness_km
from planet_reviews.crust_review import run_params
from planet_grid.sphere_graph import CubeSphere
G = CubeSphere(217)
arm = sys.argv[1]
D, TT, Z, A, L0, L1, T0, T1 = [], [], [], [], [], [], [], []
for s in range(1001, 1007):
    p = run_params(f'{arm}/seed{s}'); pl = {**p, 'lith_frame': 'lamb'}
    q = np.load(f'{arm}/seed{s}/snapshots/t+0000.0.npz')
    c = q['continental'].astype(bool); la = q['lithosphere_age_myr'].astype(float)
    m = c & np.isfinite(la)
    t0 = elevation_term(la[m], p); t1 = elevation_term(la[m], pl)
    z = q['elevation_m'][m]
    wet = z < 0
    f = np.where(wet, p['rho_mantle']/(p['rho_mantle']-p['rho_water']), 1.)
    D.append((t1-t0)*f); Z.append(z); TT.append((q['crust_thickness_m'][m]+q['sediment_m'][m])/1e3)
    A.append(G.area[m])
    L0.append(thickness_km(la[m], p)); L1.append(thickness_km(la[m], pl)); T0.append(t0); T1.append(t1)
D, Z, TT, A, L0, L1, T0, T1 = map(np.concatenate, (D, Z, TT, A, L0, L1, T0, T1))
print(f'{arm}: six seeds, 0 Myr, {D.size} continental nodes')
print(f'  C term area mean on the continent: own frame {(T0*A).sum()/A.sum():.0f} m, lamb {(T1*A).sum()/A.sum():.0f} m')
print('  crust bin | share | L own / lamb (mean km) | term change p50 / mean (m) | z p50 now / estimate')
for a, b in [(20, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]:
    k = (TT >= a) & (TT < b)
    if k.sum() > 50:
        print(f'  {a}-{b} | {A[k].sum()/A.sum():.3f} | {L0[k].mean():.0f} / {L1[k].mean():.0f} | '
              f'{np.median(D[k]):.0f} / {D[k].mean():.0f} | {np.median(Z[k]):.0f} / {np.median(Z[k]+D[k]):.0f}')
land = Z > 0; land2 = Z+D > 0
print(f'  land (continental nodes above sea) p25 / p50 / p90: now {np.round(np.percentile(Z[land], [25, 50, 90]))}'
      f', estimate {np.round(np.percentile((Z+D)[land2], [25, 50, 90]))}')
