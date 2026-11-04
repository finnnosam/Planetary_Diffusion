"""S, the margin wedge (sixty-first pass, journal 83): where the ocean-floor sediment lies, how much is subducted per
km of trench, and the outcome readings, for each arm of a sweep. Reads saved runs only; changes nothing.

  python -m planet_reviews.margin_wedge outputs\\planet-calibration\\sweep-20261104-000000 --json reviews\\margin_wedge.json

Tests (Earth on the n = 217 graph, pass 60 / journal 82; upper bounds where Earth's number includes pelagic sediment,
which the planet does not make):
  S1 ocean-floor sediment (bulk, label-free ocean: crust + sediment < 20 km) <= 192 Mkm^3 (GlobSed, CRUST1.0 where
     GlobSed has no data);
  S2 its share beyond 1,000 km from the label-free continent <= 0.14 (Earth 26 / 192);
  S3 mean sediment within 150 km of a trench (planet: subduction boundary nodes; Earth: Bird SUB) <= 0.75 km
     (Earth 9 Mkm^3 on 12 Mkm^2);
  S4 sediment subducted over the last 40 Myr per km of contour subduction (planet_reviews.crust_budget's reading)
     <= 32 km^3/Myr per km (Clift et al. 2009's 1.65 km^3/yr over Bird's 51.5k km of SUB).
Outcome readings (reported against the control): continental volume change (the run's ledger), label-free
continental area share (Earth 0.384), label-free volume (Earth 7,158 Mkm^3), crust < 30 km on the continent (Earth
0.190), sea level on the datum, deep floor, and the floor by distance seaward of passive edges (Earth's profile is the
wedge's own ceiling, so it is a build check there, not a test). Also the volumes the run booked on the wedge and in
the deep-sea fill (stats erosion_totals_km3), whole run and last 40 Myr.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_budget import contour_lengths, ocean_tests, seed_dirs, series_rates, snap
from planet_reviews.crust_review import column_thickness_m, run_params, sediment_bulk_m

SUBDUCTION = 2
EARTH = {'S1_ocean_sediment_Mkm3': 192., 'S2_far_share': 26./192., 'S3_trench_mean_km': .75, 'S4_subducted_per_km': 32.,
         'area_share': .384, 'volume_Mkm3': 7158., 'thin_share': .190}
BANDS = {'S1_ocean_sediment_Mkm3': 192., 'S2_far_share': .14, 'S3_trench_mean_km': .75, 'S4_subducted_per_km': 32.}
# Earth's floor seaward of passive edges (ETOPO p50, m; journal 82 section 4), bins in km
PROFILE_BINS = ((0, 50), (50, 100), (100, 150), (150, 200), (200, 300), (300, 500), (500, 700), (700, 1000))
EARTH_PROFILE = (-2083., -2827., -3302., -3582., -3948., -4412., -4633., -4761.)


def wq50(v, w):
    o = np.argsort(v)
    c = np.cumsum(w[o]); c /= c[-1]
    return float(v[o][min(np.searchsorted(c, .5), len(v)-1)])


def run_one(seed, g):
    p = run_params(seed)
    a = g.area
    stats = json.loads((Path(seed)/'stats.json').read_text())
    z = snap(seed, 0.)
    T = column_thickness_m(z, p)
    cont = T >= 20e3
    sed = sediment_bulk_m(z)
    bt = z['boundary_type']
    oc = ~cont
    dc = gl.distance_to(cont, g.nbr, g.dist, 1.1e8)/1e3
    dt = gl.distance_to(bt == SUBDUCTION, g.nbr, g.dist, 3e6)/1e3
    far_b = gl.distance_to(bt != 0, g.nbr, g.dist, 5e6)/1e3 > 300
    V = lambda m: float((sed[m]*a[m]).sum()/1e15)
    s1 = V(oc)
    near_t = oc & (dt < 150)
    out = {'S1_ocean_sediment_Mkm3': s1, 'S2_far_share': V(oc & (dc >= 1000))/max(s1, 1e-9),
           'S3_trench_mean_km': V(near_t)/max(a[near_t].sum()/1e12, 1e-9),
           'ocean_sediment_by_distance_Mkm3': [V(oc & (dc >= lo) & (dc < hi)) for lo, hi in
                                               ((0, 200), (200, 500), (500, 1000), (1000, 1e9))]}
    T0 = stats['series'][0]['time_myr']
    rec = series_rates(stats, -40., 0.)
    whole = series_rates(stats, T0, 0.)
    km = float(np.mean([contour_lengths(seed, float(t), g).get('subduction_km', 0.) for t in (-40, -30, -20, -10, 0)]))
    sub = rec.get('erosion_totals_km3.sediment_subducted_m3', 0.)
    out['S4_subducted_per_km'] = sub*1e6/max(km, 1.)
    out['sediment_subducted_km3_yr'] = sub
    out['subduction_contour_km'] = km
    for name, r in (('whole', whole), ('last40', rec)):
        out[f'wedge_km3_yr_{name}'] = r.get('erosion_totals_km3.margin_wedge_m3', 0.)
        out[f'deep_fill_km3_yr_{name}'] = r.get('erosion_totals_km3.deep_sea_sediment_m3', 0.)
        out[f'from_continents_km3_yr_{name}'] = r.get('erosion_totals_km3.ocean_sediment_m3', 0.)
    out['volume_change'] = stats['continental_volume_change']
    out['area_share'] = float(a[cont].sum()/a.sum())
    out['volume_Mkm3'] = float((T[cont]*a[cont]).sum()/1e15)
    out['thin_share'] = float(a[cont & (T < 30e3)].sum()/max(a[cont].sum(), 1.))
    ot = ocean_tests(z, g, p)
    out['sea_level_on_datum_m'] = ot['sea_level_on_datum_m']
    out['floor_mean_m'] = ot['floor_mean_m']
    zs = z['elevation_m'].astype(np.float64)
    out['passive_floor_p50_m'] = [wq50(zs[m], a[m]) if m.any() else float('nan') for m in
                                  (oc & far_b & (dc > lo) & (dc <= hi) for lo, hi in PROFILE_BINS)]
    return out


def mean(rs):
    m = {}
    for k in rs[0]:
        v = [r[k] for r in rs]
        m[k] = np.nanmean(np.array(v, float), 0).tolist() if isinstance(v[0], list) else float(np.mean(v))
    return m


def report(name, m, rs):
    ok = {k: m[k] <= b for k, b in BANDS.items()}
    f = lambda k: ' '.join('%.3g' % r[k] for r in rs)
    print(f'\n== {name} ({len(rs)} seeds)')
    for k in BANDS:
        print(f"  {k}: {m[k]:.3g} (<= {BANDS[k]:.3g}; Earth {EARTH[k]:.3g}) {'pass' if ok[k] else 'FAIL'}   [{f(k)}]")
    print('  ocean sediment by distance from the continent 0-200 / 200-500 / 500-1000 / > 1000 km: '
          + ' / '.join('%.0f' % x for x in m['ocean_sediment_by_distance_Mkm3']) + ' Mkm3 (Earth 85 / 51 / 29 / 26)')
    print(f"  sediment subducted, last 40 Myr: {m['sediment_subducted_km3_yr']:.2f} km3/yr over {m['subduction_contour_km']/1e3:.1f}k km")
    print(f"  wedge / deep fill / from the continents, km3/yr: whole run {m['wedge_km3_yr_whole']:.2f} / "
          f"{m['deep_fill_km3_yr_whole']:.2f} / {m['from_continents_km3_yr_whole']:.2f}; last 40 Myr "
          f"{m['wedge_km3_yr_last40']:.2f} / {m['deep_fill_km3_yr_last40']:.2f} / {m['from_continents_km3_yr_last40']:.2f}")
    print(f"  volume change {m['volume_change']*100:+.1f}%   [{f('volume_change')}]; label-free area {m['area_share']:.3f} "
          f"(Earth 0.384); volume {m['volume_Mkm3']:.0f} Mkm3 (Earth 7158); crust < 30 km {m['thin_share']:.3f} (Earth 0.190)")
    print(f"  sea level on the datum {m['sea_level_on_datum_m']:.0f} m; deep floor {m['floor_mean_m']:.0f} m")
    print('  floor seaward of passive edges, p50 m by km ' + ' '.join('%d-%d' % b for b in PROFILE_BINS))
    print('    planet ' + ' '.join('%6.0f' % x for x in m['passive_floor_p50_m']))
    print('    Earth  ' + ' '.join('%6.0f' % x for x in EARTH_PROFILE))
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='+', help='a sweep folder, combination folders, or seed folders')
    ap.add_argument('--json', help='write the results here')
    a = ap.parse_args(argv)
    g = CubeSphere(217)
    out = {'earth': EARTH, 'bands': BANDS}
    for name, seeds in seed_dirs(a.dirs):
        rs = [run_one(s, g) for s in seeds]
        m = mean(rs)
        out[name] = {'mean': m, 'passes': report(name, m, rs), 'seeds': {Path(s).name: r for s, r in zip(seeds, rs)}}
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, default=float))
    return out


if __name__ == '__main__':
    main()
