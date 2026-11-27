"""Geography of the land against Earth on the same grid: the target measures of STATE (moved here from
outputs/planet-calibration/precheck-p107-transfer/geography.py on 26 Nov, journal 149; the definitions are pass 107's,
journal 144, from the shoreline trace's class Geo).

  python -m planet_reviews.geography <arm folder> [<arm folder> ...] [--window 40] [--json out.json]
  python -m planet_reviews.geography <sweep folder> --control <arm name> [--json out.json]

Given a sweep folder, every arm in it is read (the arms are its sub-folders holding seed* folders), the control first.

Per arm and seed: the 0 Myr snapshot and the mean over the saved snapshots -window ... 0 Myr. Measures (eight-neighbour
graph of the n = 217 cube; elevation relative to sea level):
  land          land share of the globe
  holes         one-cell water holes (a water node whose eight neighbours are land), per globe
  islands       one-cell islands, per globe
  crossings     unique land-water edges (the graph's axis edges) per Mkm^2
  dry4          land whose four-ring ball is all land, % of land
  lowflat4      that dry land with centre < 500 m and four-ring relief < 300 m, % of all land
  low200        land below 200 m, % of land
  far300        land >= 300 km from water (graph distance to the nearest water node next to land), % of land
  low_far300    land below 200 m among that far land, %
The first arm is the reference of the paired differences (per seed, then the mean and sd over seeds). Earth: ETOPO on
the same grid (outputs/planet-reference/earth-n217.npz) with the ice regions (planet_earth.reference.ice_mask) left out
of every count and area (dry4 also needs its ball outside them). Earth's far-land numbers use the graph distance like
the planet's (the review's 3.8 / 19.2% used another distance and are not comparable).
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_grid.paths import OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl

KEYS = ('land', 'holes', 'islands', 'crossings', 'dry4', 'lowflat4', 'low200', 'far300', 'low_far300')


class Geo:
    def __init__(self, g, valid=None):
        self.g = g
        nbr = g.nbr.copy()
        bad = nbr < 0
        nbr[bad] = np.broadcast_to(np.arange(g.size)[:, None], nbr.shape)[bad]
        self.nb = nbr
        i, k = np.where(g.axis & (np.arange(g.size)[:, None] < np.where(bad, -1, g.nbr)))
        self.ei, self.ej = i, g.nbr[i, k]
        self.valid = np.ones(g.size, bool) if valid is None else valid
        self.A = g.area[self.valid].sum()

    def measures(self, z):
        g, a, nb, v = self.g, self.g.area, self.nb, self.valid
        land = z > 0
        hole = ~land & land[nb].all(1) & v
        island = land & (~land[nb]).all(1) & v
        L = land & v
        aland = a[L].sum()
        low, high, ok = z.copy(), z.copy(), v.copy()
        for _ in range(4):
            low = low[nb].min(1); high = high[nb].max(1); ok = ok[nb].all(1)
        dry4 = L & ok & (low > 0)
        lowflat = dry4 & (high-low < 300.) & (z < 500.)
        coast_w = ~land & land[nb].any(1)
        d = gl.distance_to(coast_w, g.nbr, g.dist, 301e3)
        far = L & (d >= 300e3)
        ev = v[self.ei] & v[self.ej]
        return {'land': float(aland/self.A), 'holes': int(hole.sum()), 'islands': int(island.sum()),
                'crossings': float(((land[self.ei] != land[self.ej]) & ev).sum()/(self.A/1e12)),
                'dry4': 100*float(a[dry4].sum()/aland), 'lowflat4': 100*float(a[lowflat].sum()/aland),
                'low200': 100*float(a[L & (z < 200.)].sum()/aland), 'far300': 100*float(a[far].sum()/aland),
                'low_far300': 100*float(a[far & (z < 200.)].sum()/max(a[far].sum(), 1.))}


def seed_read(geo, sd, window):
    snaps = sorted((sd/'snapshots').glob('t*.npz'))
    tags = {}
    for f in snaps:
        name = f.name
        if name.endswith('-plates.npz'):
            continue
        t = float(name[1:].split('.npz')[0])
        if -window-1e-6 <= t <= 1e-6:
            tags[t] = f
    if 0.0 not in tags:
        return None
    end = geo.measures(np.load(tags[0.0])['elevation_m'].astype(np.float64))
    win = [geo.measures(np.load(f)['elevation_m'].astype(np.float64)) for _, f in sorted(tags.items())]
    wmean = {k: float(np.mean([w[k] for w in win])) for k in KEYS}
    return {'end': end, 'window': wmean, 'n_window': len(win)}


def fmt(v, k):
    return f'{v:.0f}' if k in ('holes', 'islands') else (f'{v:.4f}' if k == 'land' else f'{v:.2f}')


def expand_arms(paths, control=None):
    """Arm folders from the paths given: an arm (it holds seed* folders) stays; a sweep folder becomes its arms, the
    control first."""
    arms = []
    for p in map(Path, paths):
        if any(p.glob('seed*')):
            arms.append(p)
            continue
        sub = sorted(d for d in p.iterdir() if d.is_dir() and any(d.glob('seed*')))
        if control is not None:
            sub.sort(key=lambda d: d.name != control)
        arms.extend(sub)
    if not arms:
        raise SystemExit(f'no arm folders (with seed* inside) in {paths}')
    return [str(x) for x in arms]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('arms', nargs='+')
    ap.add_argument('--window', type=float, default=40.)
    ap.add_argument('--json', default=None)
    ap.add_argument('--control', default=None, help='with a sweep folder: the arm read first')
    a = ap.parse_args(argv)
    a.arms = expand_arms(a.arms, a.control)
    g = CubeSphere(217)
    geo = Geo(g)
    out = {}
    earth_f = OUTPUTS/'planet-reference'/'earth-n217.npz'
    if earth_f.exists():
        from planet_earth.reference import ice_mask
        out['earth'] = Geo(g, ~ice_mask(g)).measures(np.load(earth_f)['elevation_m'].astype(np.float64))
    for arm in a.arms:
        arm = Path(arm)
        out[arm.name] = {sd.name: r for sd in sorted(arm.glob('seed*')) if (r := seed_read(geo, sd, a.window))}
    names = [Path(x).name for x in a.arms]
    print(f'Geography (STATE 4a; trace definitions). Window: saved snapshots -{a.window:.0f} ... 0 Myr.')
    print('measure | ' + ' | '.join(f'{n}: 0 Myr; window' for n in names) + ' | Earth')
    for k in KEYS:
        row = []
        for n in names:
            S = out[n].values()
            row.append(f"{fmt(np.mean([s['end'][k] for s in S]), k)}; {fmt(np.mean([s['window'][k] for s in S]), k)}")
        e = fmt(out['earth'][k], k) if 'earth' in out else '-'
        print(f'{k} | ' + ' | '.join(row) + f' | {e}')
    ref = names[0]
    for n in names[1:]:
        common = sorted(set(out[ref]) & set(out[n]))
        print(f'\npaired {n} minus {ref} over seeds {", ".join(s[4:] for s in common)}: mean (sd) at 0 Myr; window')
        for k in KEYS:
            de = [out[n][s]['end'][k]-out[ref][s]['end'][k] for s in common]
            dw = [out[n][s]['window'][k]-out[ref][s]['window'][k] for s in common]
            print(f'{k} | {np.mean(de):+.4g} ({np.std(de):.3g}); {np.mean(dw):+.4g} ({np.std(dw):.3g}) | per seed window: '
                  + ' '.join(f'{x:+.3g}' for x in dw))
    for n in names:
        print(f'\nper seed, {n} (0 Myr; window):')
        for s, r in out[n].items():
            print(f'{s} | ' + ' | '.join(f"{k} {fmt(r['end'][k], k)}; {fmt(r['window'][k], k)}" for k in KEYS))
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
