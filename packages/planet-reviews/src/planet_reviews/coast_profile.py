"""Land elevation by distance inland, planet against Earth (fifty-third / fifty-fourth passes).

  python -m planet_reviews.coast_profile outputs\\planet-calibration\\sweep-20261101-120000
  python -m planet_reviews.coast_profile --latest

Read-only. For every combination folder (or seed folder): at present day, on the plate grid (cube n = 217), the
great-circle distance of each land node to the nearest ocean node, and in the bins BINS_KM the area-weighted median
land elevation, the share of land below 200 m, and the median crust + sediment thickness (label-free continent,
crust_review). Seed means. Earth: ETOPO on the same grid (planet_earth.reference, ice sheets masked) and CRUST1.0
crsthk (nearest 1 degree cell).

Two numbers summarise the shape (journal 75): the rise, median land elevation at 600-1200 km inland minus at
0-50 km (Earth 387 m on n = 217), and the share of land below 200 m beyond 300 km inland (Earth 0.20). A domed
continent rises far more than Earth and has no interior lowland. Test fixed 1 Nov (pass 54): rise <= 2 x Earth's.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import DATA, OUTPUTS
from planet_reviews.crust_review import run_params, total_thickness_m, label_free_continent

BINS_KM = (0., 50., 150., 300., 600., 1200., 1e9)
RISE_BINS = (0, 4)              # bin indices: 0-50 km and 600-1200 km
RISE_FACTOR_MAX = 2.            # the test: planet rise <= 2 x Earth's
LOW_M, FAR_KM = 200., 300.


def distance_to_ocean_km(g, land):
    """Great-circle distance (km) from every land node to the nearest ocean node next to land (0 on ocean)."""
    nb = g.nbr; ok = nb >= 0
    coast = (~land) & (ok & land[np.where(ok, nb, 0)]).any(1)
    d = np.zeros(g.size)
    if not coast.any():
        d[land] = np.inf
        return d
    oc = g.xyz[coast]; idx = np.flatnonzero(land)
    for i in range(0, idx.size, 2000):
        c = g.xyz[idx[i:i+2000]] @ oc.T
        d[idx[i:i+2000]] = np.arccos(np.clip(c.max(1), -1., 1.))*g.radius/1e3
    return d


def _wmed(v, w):
    if v.size == 0:
        return np.nan
    o = np.argsort(v); c = np.cumsum(w[o])
    return float(v[o][np.searchsorted(c, c[-1]/2)])


def profile(g, z, land, T=None, cont=None, keep=None):
    """Per bin: land area share, median elevation, share below LOW_M, median thickness (km); and the rise and the
    lowland share beyond FAR_KM."""
    keep = np.ones(g.size, bool) if keep is None else keep
    land = land & keep
    d = distance_to_ocean_km(g, land)
    a = g.area; tot = a[land].sum()
    rows = []
    for lo, hi in zip(BINS_KM[:-1], BINS_KM[1:]):
        m = land & (d >= lo) & (d < hi)
        r = {'bin_km': [lo, hi], 'area_share': float(a[m].sum()/tot), 'z_p50': _wmed(z[m], a[m]),
             'lt200': float(a[m & (z < LOW_M)].sum()/max(a[m].sum(), 1.))}
        if T is not None:
            mt = m & cont
            r['T_p50_km'] = _wmed(T[mt], a[mt])/1e3
        rows.append(r)
    far = land & (d >= FAR_KM)
    return {'bins': rows, 'rise_m': rows[RISE_BINS[1]]['z_p50']-rows[RISE_BINS[0]]['z_p50'],
            'lowland_far': float(a[far & (z < LOW_M)].sum()/max(a[far].sum(), 1.))}


def earth(g):
    from planet_earth.reference import ice_mask
    z = np.load(OUTPUTS/'planet-reference'/f'earth-n{g.n}.npz')['elevation_m'].astype(float)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)*1e3
    lat = np.degrees(np.arcsin(g.xyz[:, 2])); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    Te = T[np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)]
    return profile(g, z, z > 0, Te, label_free_continent(Te), ~ice_mask(g))


def planet(g, seed_dir):
    z = np.load(Path(seed_dir)/'snapshots'/'t+0000.0.npz')
    e = z['elevation_m'].astype(float)          # relative to sea level already (sea_level_m is the sea on the datum)
    T = total_thickness_m(z, run_params(seed_dir))
    return profile(g, e, e > 0, T, label_free_continent(T))


def mean_profile(ps):
    out = {'bins': [], 'rise_m': float(np.mean([p['rise_m'] for p in ps])),
           'lowland_far': float(np.mean([p['lowland_far'] for p in ps]))}
    for k in range(len(BINS_KM)-1):
        rs = [p['bins'][k] for p in ps]
        out['bins'].append({key: (rs[0][key] if key == 'bin_km' else float(np.nanmean([r[key] for r in rs])))
                            for key in rs[0]})
    return out


def show(name, p, e=None):
    print(f"{name}: rise {p['rise_m']:.0f} m, land < {LOW_M:.0f} m beyond {FAR_KM:.0f} km {p['lowland_far']:.3f}"
          + (f" (Earth {e['rise_m']:.0f} m, {e['lowland_far']:.3f}; test rise <= {RISE_FACTOR_MAX:.0f} x Earth: "
             f"{'pass' if p['rise_m'] <= RISE_FACTOR_MAX*e['rise_m'] else 'fail'})" if e else ''))
    for r in p['bins']:
        lo, hi = r['bin_km']
        print(f"    {lo:5.0f}-{hi if hi < 1e8 else float('inf'):<6.0f} km: area {r['area_share']:.3f}, z p50 {r['z_p50']:6.0f} m, "
              f"< {LOW_M:.0f} m {r['lt200']:.3f}" + (f", crust p50 {r['T_p50_km']:.1f} km" if 'T_p50_km' in r else ''))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='*', help='a sweep folder, combination folders, or seed folders')
    ap.add_argument('--latest', action='store_true', help='the newest outputs/planet-calibration/sweep-* folder')
    ap.add_argument('--json', default=None)
    a = ap.parse_args(argv)
    if a.latest:
        a.dirs = [str(sorted((OUTPUTS/'planet-calibration').glob('sweep-*'))[-1])]
        print('sweep', a.dirs[0])
    g = CubeSphere(217)
    e = earth(g)
    show('Earth (ETOPO, no ice sheets; CRUST1.0)', e)
    combos = []
    for d in map(Path, a.dirs):
        if (d/'snapshots').is_dir():
            combos.append((d.parent.name+'/'+d.name, [d]))
        elif any(d.glob('seed*/snapshots')):
            combos.append((d.name, sorted(d.glob('seed*'))))
        else:
            combos += [(c.name, sorted(c.glob('seed*'))) for c in sorted(d.iterdir()) if c.is_dir() and any(c.glob('seed*/snapshots'))]
    out = {'earth': e}
    for name, seeds in combos:
        ps = [planet(g, s) for s in seeds]
        m = mean_profile(ps)
        m['seeds'] = {s.name: {'rise_m': p['rise_m'], 'lowland_far': p['lowland_far']} for s, p in zip(seeds, ps)}
        m['passes'] = bool(m['rise_m'] <= RISE_FACTOR_MAX*e['rise_m'])
        out[name] = m
        show(f'{name} ({len(ps)} seeds)', m, e)
        print('    by seed: ' + ', '.join(f"{k} {v['rise_m']:.0f} m / {v['lowland_far']:.3f}" for k, v in m['seeds'].items()))
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, default=float))
    return out


if __name__ == '__main__':
    main()
