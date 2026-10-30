"""Field-by-field comparison of two runs' snapshots (identity checks between sweeps; NaN-aware).

  python tools\\compare_runs.py <seed folder A> <seed folder B>
  python tools\\compare_runs.py <combination folder A> <combination folder B>   (every seedNNNN present in both)
  python tools\\compare_runs.py A B --ignore extension_age_myr,thickening_age_myr      (skip fields only one side has)

Prints, per seed, the snapshots compared, the fields that differ, and the stats.json keys that differ.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def compare(a, b, ignore=()):
    sa = sorted((a/'snapshots').glob('t*.npz')); sb = {x.name: x for x in (b/'snapshots').glob('t*.npz')}
    diff, n = [], 0
    for x in sa:
        if x.name not in sb:
            diff.append((x.name, 'missing in B')); continue
        za, zb = np.load(x), np.load(sb[x.name])
        for k in sorted((set(za.files) | set(zb.files)) - set(ignore)):
            n += 1
            if k not in za.files or k not in zb.files:
                diff.append((x.name, k+' (only one side)')); continue
            u, v = za[k], zb[k]
            if u.shape != v.shape or not np.array_equal(u, v, equal_nan=u.dtype.kind == 'f'):
                diff.append((x.name, k))
    ja = json.loads((a/'stats.json').read_text()); jb = json.loads((b/'stats.json').read_text())
    keys = sorted(k for k in set(ja) | set(jb) if ja.get(k) != jb.get(k))
    print(f"{a.name}: {len(sa)} snapshots, {n} fields{f' (ignoring {sorted(ignore)})' if ignore else ''}, {len(diff)} differ{': '+', '.join(f'{t} {k}' for t, k in diff[:8]) if diff else ''}; "
          f"stats.json keys that differ: {keys}")
    return not diff


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('--ignore', default='', help='comma-separated snapshot fields to leave out of the comparison')
    o = ap.parse_args(argv)
    A, B, ignore = Path(o.a), Path(o.b), tuple(f for f in o.ignore.split(',') if f)
    pairs = [(A, B)] if (A/'snapshots').is_dir() else [(A/s.name, B/s.name) for s in sorted(A.glob('seed*')) if (B/s.name).is_dir()]
    ok = all([compare(a, b, ignore) for a, b in pairs])
    print('IDENTICAL' if ok else 'DIFFERENT', f'({len(pairs)} seeds)')


if __name__ == '__main__':
    main()
