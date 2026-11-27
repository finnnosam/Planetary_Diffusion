"""Pass 55: the planet's height at total thickness, split by provenance (interior_probe's definition: collided =
orogeny clock reset during the run; inherited = neither collided nor new) and, for collided crust, by orogeny age.
Relative to sea level and on the datum (z + sea_level_m). Six seeds, 0 Myr; area shares of the label-free continent."""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import run_params, total_thickness_m
BINS = [30, 33, 36, 39, 42, 45, 50, 60]
def wmed(v, w):
    if len(v) < 5: return np.nan
    o = np.argsort(v); c = np.cumsum(w[o]); return float(v[o][np.searchsorted(c, c[-1]/2)])
g = CubeSphere(217); a = g.area
sw = Path(sys.argv[1])
for combo in sys.argv[2:]:
    acc = {}
    for sd in sorted((sw/combo).glob('seed*')):
        z = np.load(sd/'snapshots'/'t+0000.0.npz'); p = run_params(sd)
        T = total_thickness_m(z, p)/1e3; cont = T >= 20
        e = z['elevation_m'].astype(float) + float(z['sea_level_m'])          # on the datum
        oa = z['orogeny_age_myr'].astype(float); new = z['new_crust'].astype(bool)
        col = np.isfinite(oa) & (oa < 250)
        groups = {'inherited': cont & ~col & ~new, 'collided': cont & col & ~new,
                  'coll <30Myr': cont & col & ~new & (oa < 30), 'coll >=30Myr': cont & col & ~new & (oa >= 30)}
        for k, m in groups.items():
            r = acc.setdefault(k, {'z': [], 'share': []})
            r['z'].append([wmed(e[m & (T >= lo) & (T < hi)], a[m & (T >= lo) & (T < hi)]) for lo, hi in zip(BINS[:-1], BINS[1:])])
            r['share'].append([a[m & (T >= lo) & (T < hi)].sum()/a[cont].sum() for lo, hi in zip(BINS[:-1], BINS[1:])])
    print(f'{combo}: z on the datum (m) by total thickness; area share of the continent in brackets')
    print('   total km      ', ' '.join(f'{lo:>2}-{hi:<9}' for lo, hi in zip(BINS[:-1], BINS[1:])))
    for k, r in acc.items():
        zz = np.nanmean(np.array(r['z'], float), 0); ss = np.mean(r['share'], 0)
        print(f'   {k:<14}', ' '.join(f'{v:5.0f} ({s:.3f})' for v, s in zip(zz, ss)))
