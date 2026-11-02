"""Land by distance inland, against the Airy line of its own column (fifty-fifth pass's probe, promoted in the
fifty-sixth; test C1 of C, the continental lithosphere's thermal state).

In coast_profile's distance bins (BINS_KM), per band of land: the area-weighted median elevation on the datum (Earth:
ETOPO, sea level = datum; planet: elevation_m + sea_level_m), crystalline crust, sediment, the planet's Airy line of
that column (cont_ref_thick_m at cont_ref_elev_m, rho_crust against rho_mantle; sediment at rho_sediment, bulk) and
z - Airy. Distance bins do not depend on CRUST1.0's thickness errors, so band medians are free of the regression
dilution that biases bins by thickness (pass 34). Earth: CRUST1.0 crsthk and sedthk (nearest 1-degree cell), ice
sheets excluded.

C1 (fixed 1 Nov, journal 77 section 6): the fall of z - Airy from the coast band (0-50 km) to the interior band
(600-1200 km). Earth -773 m on n = 217: margins stand ~400 m above the Airy line and old interiors ~370 m below it
(old, cold mantle lithosphere). Band: between 0.5 and 1.5 x Earth's (-390 ... -1,160 m). The planet's lithosphere age
(lithosphere_age_myr, C) is printed by band where the snapshot has it; Earth's from CRUST1.0's types
(planet_tectonics.lithosphere.CLASS_AGE_MA).

Reported beside C1, not scored: the planet's z - Airy without its dynamic topography (z - dynamic_topography_m), C1
read on it ('c1_nodyn_m'). The planet's random field (fitted dyn_std_m 500, power to degree 30) puts interior land on
its highs (+150-240 m at 300-1200 km inland on the pass-56 pre-check, ~0 at the coast), where Earth's field has 86% of
its power at degrees 1-3 (Valentine & Davies 2020) and cannot make a coast-to-interior contrast of ~600 km.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.coast_profile import BINS_KM, RISE_BINS, distance_to_ocean_km
from planet_reviews.crust_review import run_params, sediment_bulk_m

C1_FACTORS = (.5, 1.5)             # the band: 0.5 ... 1.5 x Earth's fall


def wmed(v, w):
    ok = np.isfinite(v)
    v, w = v[ok], w[ok]
    if w.sum() == 0:
        return np.nan
    o = np.argsort(v); c = np.cumsum(w[o])
    return float(v[o][np.searchsorted(c, c[-1]/2)])


def airy_m(Tc, S, p):
    """The planet's Airy line (m on the datum) of a column of crystalline crust Tc and bulk sediment S (m)."""
    rm = p['rho_mantle']
    return p['cont_ref_elev_m']+(Tc-p['cont_ref_thick_m'])*(1-p['rho_crust']/rm)+S*(1-p['rho_sediment']/rm)


def budget(g, z, Tc, S, land, p, age=None, dyn=None):
    """Band medians by distance inland: z, crystalline km, sediment, Airy, z - Airy (m), lithosphere age (Myr), and
    the coast-to-interior changes (rise of z, of the Airy line, and C1, the fall of z - Airy; with dyn, also z - dyn -
    Airy and its C1)."""
    a = g.area
    d = distance_to_ocean_km(g, land)
    ai = airy_m(Tc, S, p)
    rows = []
    for lo, hi in zip(BINS_KM[:-1], BINS_KM[1:]):
        m = land & (d >= lo) & (d < hi)
        rows.append({'bin_km': [lo, hi], 'area_share': float(a[m].sum()/max(a[land].sum(), 1.)),
                     'z': wmed(z[m], a[m]), 'crystalline_km': wmed(Tc[m], a[m])/1e3, 'sediment': wmed(S[m], a[m]),
                     'airy': wmed(ai[m], a[m]), 'resid': wmed((z-ai)[m], a[m]),
                     'age_myr': wmed(age[m], a[m]) if age is not None else np.nan,
                     'resid_nodyn': wmed((z-dyn-ai)[m], a[m]) if dyn is not None else np.nan})
    i, j = RISE_BINS
    return {'bins': rows, 'rise_m': rows[j]['z']-rows[i]['z'], 'airy_rise_m': rows[j]['airy']-rows[i]['airy'],
            'c1_m': rows[j]['resid']-rows[i]['resid'], 'c1_nodyn_m': rows[j]['resid_nodyn']-rows[i]['resid_nodyn']}


def earth(g, p):
    from planet_earth.reference import ice_mask
    from planet_tectonics.lithosphere import CLASS_AGE_MA
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    tn = np.array([l.split() for l in (DATA/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])
    age = np.array([CLASS_AGE_MA.get(t, np.nan) for t in tn[np.clip((89.5-lat+.5).astype(int), 0, 179), jj]])
    z = np.load(OUTPUTS/'planet-reference'/f'earth-n{g.n}.npz')['elevation_m'].astype(float)
    land = (z > 0) & ~ice_mask(g)
    return budget(g, z, T-S, S, land, p, age)


def planet(g, seed_dir, t='t+0000.0'):
    z = np.load(Path(seed_dir)/'snapshots'/f'{t}.npz')
    p = run_params(seed_dir)
    e = z['elevation_m'].astype(float)             # relative to sea level (sea_level_m is the sea on the datum)
    c = z['continental'].astype(bool)
    S = np.where(c, sediment_bulk_m(z), 0.)
    age = z['lithosphere_age_myr'].astype(float) if 'lithosphere_age_myr' in z.files else None
    dyn = z['dynamic_topography_m'].astype(float) if 'dynamic_topography_m' in z.files else None
    return budget(g, e+float(z['sea_level_m']), z['crust_thickness_m'].astype(float), S, e > 0, p, age, dyn)


def mean_budget(bs):
    out = {k: float(np.nanmean([b[k] for b in bs])) for k in ('rise_m', 'airy_rise_m', 'c1_m', 'c1_nodyn_m')}
    out['bins'] = [{key: (bs[0]['bins'][k][key] if key == 'bin_km' else float(np.nanmean([b['bins'][k][key] for b in bs])))
                    for key in bs[0]['bins'][k]} for k in range(len(BINS_KM)-1)]
    return out


def c1_band(e):
    return sorted((C1_FACTORS[0]*e['c1_m'], C1_FACTORS[1]*e['c1_m']))


def show(name, b, e=None):
    lo, hi = c1_band(e) if e else (np.nan, np.nan)
    verdict = ('pass' if lo <= b['c1_m'] <= hi else 'fail') if e else ''
    print(f"{name}: rise z {b['rise_m']:.0f} m, Airy {b['airy_rise_m']:.0f} m, C1 (z - Airy, interior minus coast) "
          f"{b['c1_m']:.0f} m" + (f"  [band {lo:.0f} ... {hi:.0f}: {verdict}]" if e else '')
          + (f"; without the dynamic topography {b['c1_nodyn_m']:.0f} m (reported)" if np.isfinite(b.get('c1_nodyn_m', np.nan)) else ''))
    for r in b['bins']:
        lo_, hi_ = r['bin_km']
        print(f"    {lo_:5.0f}-{hi_ if hi_ < 1e8 else float('inf'):<6.0f} km: z {r['z']:6.0f}  crystalline "
              f"{r['crystalline_km']:5.1f} km  sediment {r['sediment']:5.0f}  Airy {r['airy']:6.0f}  z - Airy "
              f"{r['resid']:6.0f}" + (f"  (without dyn {r['resid_nodyn']:6.0f})" if np.isfinite(r['resid_nodyn']) else '')
              + (f"  lithosphere age {r['age_myr']:6.0f} Myr" if np.isfinite(r['age_myr']) else ''))


def main(argv=None):
    from planet_tectonics.plates import PARAMS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='*', help='a sweep folder, combination folders, or seed folders')
    ap.add_argument('--latest', action='store_true', help='the newest outputs/planet-calibration/sweep-* folder')
    ap.add_argument('--json', default=None)
    a = ap.parse_args(argv)
    if a.latest:
        a.dirs = [str(sorted((OUTPUTS/'planet-calibration').glob('sweep-*'))[-1])]
        print('sweep', a.dirs[0])
    g = CubeSphere(217)
    e = earth(g, PARAMS)
    show('Earth (ETOPO, no ice sheets; CRUST1.0; the planet\'s Airy line)', e)
    combos = []
    for d in map(Path, a.dirs):
        if (d/'snapshots').is_dir():
            combos.append((d.parent.name+'/'+d.name, [d]))
        elif any(d.glob('seed*/snapshots')):          # a combination folder: named with its sweep
            combos.append((d.parent.name+'/'+d.name, sorted(d.glob('seed*'))))
        else:
            combos += [(c.name, sorted(c.glob('seed*'))) for c in sorted(d.iterdir()) if c.is_dir() and any(c.glob('seed*/snapshots'))]
    out = {'earth': e, 'c1_band_m': c1_band(e)}
    for name, seeds in combos:
        bs = [planet(g, s) for s in seeds]
        m = mean_budget(bs)
        m['seeds'] = {s.name: {'c1_m': b['c1_m'], 'c1_nodyn_m': b['c1_nodyn_m'], 'rise_m': b['rise_m']}
                      for s, b in zip(seeds, bs)}
        lo, hi = out['c1_band_m']
        m['c1_passes'] = bool(lo <= m['c1_m'] <= hi)
        out[name] = m
        show(f'{name} ({len(bs)} seeds, means of band medians)', m, e)
        print('    by seed: C1 ' + ', '.join(f"{k} {v['c1_m']:.0f}" for k, v in m['seeds'].items())
              + ' | without dyn ' + ', '.join(f"{k} {v['c1_nodyn_m']:.0f}" for k, v in m['seeds'].items()))
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, default=float))
    return out


if __name__ == '__main__':
    main()
