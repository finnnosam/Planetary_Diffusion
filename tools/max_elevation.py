"""Highest node elevation in every saved snapshot of a sweep's runs (the 'no towers' veto, journal 86).

  python tools\\max_elevation.py <sweep folder or combination folder> [--limit 8850]

Earth's highest summit is 8,849 m; on the n = 217 graph ETOPO's highest node mean is 5,747 m, so a node above the
limit is a column no Earth process holds (the no-flow arm of sweep-20261101-191306 reached 12-20 km on every seed).
Prints per combination and seed the maximum over all snapshots and when it occurred, and the verdict per combination.
"""
import argparse
from pathlib import Path

import numpy as np


def run_max(seed_dir):
    best = (-np.inf, None)
    for f in sorted(Path(seed_dir, 'snapshots').glob('t*.npz')):
        if f.name.endswith('-plates.json'):
            continue
        z = np.load(f)
        m = float(z['elevation_m'].max())
        if m > best[0]:
            best = (m, float(z['time_myr']))
    return best


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folder')
    ap.add_argument('--limit', type=float, default=8850.)
    a = ap.parse_args(argv)
    seeds = sorted(p for p in Path(a.folder).rglob('seed*') if (p/'snapshots').is_dir())
    combos = {}
    for s in seeds:
        combos.setdefault(s.parent.name, []).append(s)
    for c, ss in combos.items():
        res = [(s.name, *run_max(s)) for s in ss]
        ok = all(m <= a.limit for _, m, _ in res)
        print(f"{c}: " + '; '.join(f"{n} {m:,.0f} m at {t:+.0f} Myr" for n, m, t in res)
              + f" -> {'passes' if ok else 'FAILS'} (limit {a.limit:,.0f} m)")


if __name__ == '__main__':
    main()
