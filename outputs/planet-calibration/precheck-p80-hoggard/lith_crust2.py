"""Pass 80 (journal 104): the planet's lithosphere against its crust, each run in its own frame (lith_frame), beside
Earth. pass 77's lith_crust.py with L from lithosphere.thickness_km under the run's parameters (crust_review.run_params;
it computed Artemieva's line for every run) and Earth's node-level reference added.
Earth: Lamb et al. 2020 ds02 (11 continental means of Priestley & McKenzie 2016 seismic L on CRUST1.0 crust; Lamb's
reduced elevation z + k (L - 100), k 10.753 m per km), and Hoggard et al. 2020 SL2013sv L on the n = 217 nodes of
Earth's continent (lithosphere.earth_continent's domain): dL/dT 3.18 km per km, area means by crust bin 97 / 98 / 107 /
143 / 158 / 157 / 197 km (precheck-p80-hoggard/hoggard_law.txt), the like-for-like reading of a node-level slope.
Planet: continental nodes (crust + sediment >= 20 km) at 0 Myr, six seeds of each arm given. Run from the repo root."""
import sys, numpy as np
from planet_tectonics.lithosphere import thickness_km
from planet_reviews.crust_review import run_params
K = 10.753
med = lambda x, t, a, b: float(np.median(x[(t >= a) & (t < b)]))
rows = [l.split('\t') for l in open('data/Lamb2020/ds02.xyz') if not l.startswith('#') and l.strip()]
d = np.array([[float(x) for x in r[1:7]] for r in rows]); names = [r[0].strip() for r in rows]
z, zr, L, T = d[:, 1], d[:, 2], d[:, 3], d[:, 4]
print('Earth (Lamb et al. 2020 ds02, 11 continental means):')
for nm, s in (('all 11', np.ones(len(names), bool)), ('without Antarctica', np.array(['Antarctica' not in n for n in names]))):
    print(f'  {nm}: dL/dT {np.polyfit(T[s], L[s], 1)[0]:.1f} km/km, corr {np.corrcoef(T[s], L[s])[0, 1]:.2f}; '
          f'reduced z vs T {np.polyfit(T[s], zr[s], 1)[0]:.0f} m/km; raw z vs T {np.polyfit(T[s], z[s], 1)[0]:.0f} m/km')
print('Earth (Hoggard et al. 2020 SL2013sv, n = 217 continental nodes): dL/dT 3.18 km/km; area-mean L by crust bin '
      '20-25 97, 25-30 98, 30-35 107, 35-40 143, 40-45 158, 45-55 157, 55-80 197 km')
for arm in sys.argv[1:]:
    Z, TT, LL = [], [], []
    for s in range(1001, 1007):
        p = run_params(f'{arm}/seed{s}')
        q = np.load(f'{arm}/seed{s}/snapshots/t+0000.0.npz')
        c = q['continental'].astype(bool); la = q['lithosphere_age_myr'].astype(float)
        t = (q['crust_thickness_m']+q['sediment_m'])/1e3
        m = c & np.isfinite(la) & (t >= 20)
        Z.append(q['elevation_m'][m]); TT.append(t[m]); LL.append(thickness_km(la[m], p))
    frame = p.get('lith_frame', 'artemieva')
    Z, TT, LL = map(np.concatenate, (Z, TT, LL))
    red = Z+K*(LL-100)
    print(f'\n{arm} (frame {frame}): corr(L, T) {np.corrcoef(TT, LL)[0, 1]:.2f}, dL/dT {np.polyfit(TT, LL, 1)[0]:.2f} km/km, '
          f'reduced z vs T {np.polyfit(TT, red, 1)[0]:.0f} m/km, raw z vs T {np.polyfit(TT, Z, 1)[0]:.0f} m/km (node least '
          f'squares); bin medians 30-35 -> 40-45 km: raw {(med(Z, TT, 40, 45)-med(Z, TT, 30, 35))/10:.0f} m/km, L '
          f'{(med(LL, TT, 40, 45)-med(LL, TT, 30, 35))/10:.1f} km/km')
    print('  T bin | L p25/p50/p75 km | L mean | z p50 m | reduced z p50 m')
    for a, b in [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]:
        m = (TT >= a) & (TT < b)
        if m.sum() > 50:
            print(f'  {a}-{b} | {np.round(np.percentile(LL[m], [25, 50, 75]))} | {LL[m].mean():.0f} | '
                  f'{np.median(Z[m]):.0f} | {np.median(red[m]):.0f}')
