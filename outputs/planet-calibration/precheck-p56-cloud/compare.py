"""Pass 56 pre-check: combo-03 + C against combo-03, cloud seed 1007 (paired: same seed, same start). Read-only.
C1 (coast_budget), P1-P5 (as the Collapse Flow Sweep's rule), the guards readable in the cloud (volume, export,
plates), and reported values. 0 Myr, and the mean over the last 40 Myr's saved snapshots where noted."""
import sys, json, glob, re
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import total_thickness_m, wpct, run_params
from planet_reviews.interior_probe import slope_m_km, edge_depth_km
from planet_reviews import coast_profile as cp, coast_budget as cb

g = CubeSphere(217); A = g.area


def shape(f, p):
    z = np.load(f); e = z['elevation_m'].astype(float); land = e > 0
    T = total_thickness_m(z, p); c = T >= 20e3
    s = slope_m_km(g, e); hi = land & (e > 2000)
    thick = c & (T > 45e3)
    dep, _ = edge_depth_km(g, thick)
    pr = cp.profile(g, e, land, T, c)
    out = {'P1_gt45': A[thick].sum()/A[c].sum(), 'P2_gt55': A[c & (T > 55e3)].sum()/A[c].sum(),
           'P3_slope_hi_p50': wpct(s[hi], A[hi], 50) if hi.any() else np.nan,
           'P4_depth45_p50': wpct(dep[thick], A[thick], 50) if thick.any() else np.nan, 'P5_rise': pr['rise_m'],
           'area': A[c].sum()/A.sum(), 'land_frac': A[land].sum()/A.sum(),
           'land_p25': wpct(e[land], A[land], 25), 'land_p50': wpct(e[land], A[land], 50), 'land_p90': wpct(e[land], A[land], 90),
           'drowned': A[c & ~land].sum()/A[c].sum(), 'lowland_far': pr['lowland_far'], 'sea_datum': float(z['sea_level_m'])}
    if 'lithosphere_age_myr' in z.files:
        la = z['lithosphere_age_myr'].astype(float); cc = z['continental'].astype(bool)
        out['lith_age_p10'], out['lith_age_p50'], out['lith_age_p90'] = (wpct(la[cc], A[cc], q) for q in (10, 50, 90))
    return out


res = {}
for run in sys.argv[1:]:
    sd = glob.glob(f'{run}/seed*')[0]
    st = json.load(open(f'{sd}/stats.json')); p = run_params(sd)
    e = st['present']['erosion_totals_km3']
    export = -(e['continental_eroded_m3']+e['continental_deposited_m3']-e.get('continental_deposited_submerged_m3', 0.))/st['years']/1e6
    fs = sorted(glob.glob(f'{sd}/snapshots/t*.npz'), key=lambda f: float(re.search(r't([+-][\d.]+)\.npz', f).group(1)))
    late = [f for f in fs if float(re.search(r't([+-][\d.]+)\.npz', f).group(1)) >= -40]
    ms = [shape(f, p) for f in late]
    b = cb.planet(g, sd)
    now = ms[-1]
    res[run] = {'now': now, 'late': {k: float(np.nanmean([m[k] for m in ms])) for k in now}, 'budget': b,
                'volume': st['continental_volume_change'], 'crust_volume': st['continental_crust_volume_change'],
                'export': export, 'plates': st['present']['plates']}
    print(f'== {run}: volume {st["continental_volume_change"]:+.3f} (crust {st["continental_crust_volume_change"]:+.3f}), '
          f'export {export:.2f} km3/yr, plates at 0 Myr {st["present"]["plates"]}')
    cb.show('   coast budget, 0 Myr', b)
    for k in now:
        print(f'   {k:16s} 0 Myr {now[k]:9.3f}   last-40 mean {res[run]["late"][k]:9.3f}')
json.dump(res, open('p56_compare.json', 'w'), indent=1, default=float)
