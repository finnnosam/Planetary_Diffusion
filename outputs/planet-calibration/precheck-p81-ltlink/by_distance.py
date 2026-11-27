"""Pass 81: crust and L by distance from the continent's edge: Earth (CRUST1.0 total, SL2013sv; the lamb law on Earth's
class ages) against planet snapshots (crust + sediment; L in the run's frame). Args: snapshot dirs' arm and time.
Run from the repo root."""
import sys, numpy as np
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import lithosphere as lt
from planet_reviews.crust_review import run_params
g = CubeSphere(217); B = lt.DIST_BINS_KM
a, d, age = lt.earth_continent()
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]
tn = np.array([l.split() for l in (DATA/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
from planet_earth.reference import ice_mask
cont = ~np.isin(tn, lt.OCEANIC_TYPES) & (T >= 20) & ~ice_mask(g)
Te = T[cont]; Llaw = lt.thickness_km(age, {'lith_frame': 'lamb'})
rows = {}
def put(key, dist_km, crust, L, w):
    out = []
    for lo, hi in zip(B[:-1], B[1:]):
        k = (dist_km >= lo) & (dist_km < hi)
        out.append((w[k].sum()/w.sum(), np.average(crust[k], weights=w[k]) if k.any() else np.nan,
                    np.average(L[k], weights=w[k]) if k.any() else np.nan))
    rows[key] = out
put('Earth (law on ages)', d/1e3, Te, Llaw, a)
for arm in sys.argv[2:]:
    D, C, LL, W = [], [], [], []
    for s in range(1001, 1007):
        p = run_params(f'{arm}/seed{s}'); q = np.load(f'{arm}/seed{s}/snapshots/{sys.argv[1]}.npz')
        t = (q['crust_thickness_m']+q['sediment_m'])/1e3; la = q['lithosphere_age_myr'].astype(float)
        whole = q['continental'].astype(bool) & (t >= 20)
        dd = gl.distance_to(~whole, g.nbr, g.dist)/1e3
        m = whole & np.isfinite(la)
        D.append(dd[m]); C.append(t[m]); LL.append(lt.thickness_km(la[m], p)); W.append(g.area[m])
    put(arm.split('/')[-1]+' '+sys.argv[1], *map(np.concatenate, (D, C, LL, W)))
print('distance bin | ' + ' | '.join(f'{k}: area / crust / L' for k in rows))
for i, (lo, hi) in enumerate(zip(B[:-1], B[1:])):
    print(f'{lo:5.0f}-{hi:5.0f} | ' + ' | '.join(f'{v[i][0]:.3f} / {v[i][1]:.1f} / {v[i][2]:.0f}' for v in rows.values()))
