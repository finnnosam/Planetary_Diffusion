"""Field-by-field comparison of two runs (identity checks between sweeps; NaN-aware).

  python tools\\compare_runs.py <seed folder A> <seed folder B>
  python tools\\compare_runs.py <combination folder A> <combination folder B>   (every seedNNNN present in both)
  python tools\\compare_runs.py A B --ignore extension_age_myr,thickening_age_myr      (skip these snapshot fields)

Compares what both runs have: the snapshots saved by both (a 5 Myr sweep against a 1 Myr one shares every 5th), the
fields both snapshots carry, and every stats.json value both sides hold (recursively, through the series), except the
run's wall time ('seconds') and its command-line overrides ('params'; 'params_effective' is compared). Prints, per
seed, the value differences, and apart from them what only one side has (snapshots, fields, stats keys: new or
removed outputs and parameters). The verdict, IDENTICAL or DIFFERENT, is on the shared values only (forty-eighth pass;
before it, fields and keys present on one side only also printed DIFFERENT).
"""
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

STATS_SKIP = ('seconds', 'params')


def _walk(a, b, path, diff, only):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            p = f'{path}.{k}' if path else str(k)
            if k not in a:
                only.append((p, 'B'))
            elif k not in b:
                only.append((p, 'A'))
            else:
                _walk(a[k], b[k], p, diff, only)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (u, v) in enumerate(zip(a, b)):
            _walk(u, v, f'{path}[{i}]', diff, only)
    elif isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return
    elif a != b:
        diff.append(path)


def _collapse(paths):
    """'series[12].x' -> 'series[*].x', counted."""
    out = {}
    for p in paths:
        q = re.sub(r'\[\d+\]', '[*]', p)
        out[q] = out.get(q, 0)+1
    return [f'{k} (x{n})' if n > 1 else k for k, n in sorted(out.items())]


def compare_stats(a, b):
    ja = json.loads((a/'stats.json').read_text()); jb = json.loads((b/'stats.json').read_text())
    for k in STATS_SKIP:
        ja.pop(k, None); jb.pop(k, None)
    diff, only = [], []
    _walk(ja, jb, '', diff, only)
    return diff, only


def compare(a, b, ignore=()):
    na = {x.name for x in (a/'snapshots').glob('t*.npz')}; nb = {x.name for x in (b/'snapshots').glob('t*.npz')}
    shared = sorted(na & nb)
    diff, only, n = [], [], 0
    for name in shared:
        za, zb = np.load(a/'snapshots'/name), np.load(b/'snapshots'/name)
        for k in sorted((set(za.files) | set(zb.files))-set(ignore)):
            if k not in za.files or k not in zb.files:
                only.append((k, 'A' if k in za.files else 'B')); continue
            n += 1
            u, v = za[k], zb[k]
            if u.shape != v.shape or not np.array_equal(u, v, equal_nan=u.dtype.kind == 'f'):
                diff.append(f'{name} {k}')
    sd, so = compare_stats(a, b)
    snap_only = {'A': sorted(na-nb), 'B': sorted(nb-na)}
    field_only = {s: sorted({k for k, side in only if side == s}) for s in 'AB'}
    stats_only = {s: _collapse([p for p, side in so if side == s]) for s in 'AB'}
    ok = not diff and not sd
    print(f"{a.name}: {len(shared)} shared snapshots, {n} shared fields{f' (ignoring {sorted(ignore)})' if ignore else ''}: "
          f"{len(diff)} differ{': '+', '.join(diff[:8]) if diff else ''}; stats.json shared values that differ: "
          f"{len(sd)}{': '+', '.join(_collapse(sd)[:8]) if sd else ''}")
    for s in 'AB':
        parts = []
        if snap_only[s]:
            parts.append(f"{len(snap_only[s])} snapshots ({snap_only[s][0]} ... {snap_only[s][-1]})")
        if field_only[s]:
            parts.append('fields '+', '.join(field_only[s]))
        if stats_only[s]:
            parts.append('stats keys '+', '.join(stats_only[s][:12])+(' ...' if len(stats_only[s]) > 12 else ''))
        if parts:
            print(f"  only in {s}: "+'; '.join(parts))
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('--ignore', default='', help='comma-separated snapshot fields to leave out of the comparison')
    o = ap.parse_args(argv)
    A, B, ignore = Path(o.a), Path(o.b), tuple(f for f in o.ignore.split(',') if f)
    pairs = [(A, B)] if (A/'snapshots').is_dir() else [(A/s.name, B/s.name) for s in sorted(A.glob('seed*')) if (B/s.name).is_dir()]
    if not pairs:
        raise SystemExit(f'no seed folder in both {A} and {B}')
    ok = all([compare(a, b, ignore) for a, b in pairs])
    print('IDENTICAL' if ok else 'DIFFERENT', f'on the shared values ({len(pairs)} seeds; one-side-only items listed above)')


if __name__ == '__main__':
    main()
