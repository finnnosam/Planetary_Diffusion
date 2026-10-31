"""Probe (pass 42): where the intra-oceanic trench of the window comes from, by plate pair. For every snapshot of the
window (-100 ... 0 Myr, every 10 Myr, plus -5 ... 0), the contour's intra-oceanic and continental-side trench by
unordered plate pair, classed by the pair's origin:
  start-start      both plates from the start
  sibling:<kind>   one plate made from the other by an event (rift / subduction initiation / detach), possibly via
                   a chain of events of that plate's line
  created:<kind>   one plate created by an event (the younger plate's kind), the other not its parent
and, per pair, the first snapshot (10 Myr spacing before -40, 1 Myr after) since which it has had convergent
contour without a break ("onset").

  python trench_pair_origin.py <seed dir> [...]
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_reviews import boundary_lengths as bl
from planet_reviews import trench_origin as to

ALL_T = [float(t) for t in range(-250, -40, 10)]+[float(t) for t in range(-40, 1)]
WIN_T = [float(t) for t in range(-100, 1, 10)]


def pair_lengths(g, f):
    z = np.load(str(f)+'.npz'); plates = json.loads(Path(str(f)+'-plates.json').read_text())
    pid = z['plate_id']; cont = z['continental'].astype(bool)
    s = bl.typed_segments(g, pid, plates, cont)
    L, typ, cv, na, nb, ka, kb = s['L'], s['typ'], s['conv'], s['na'], s['nb'], s['ka'], s['kb']
    oo = ~(cont[na] | cont[nb])
    sub = typ == 2
    lo, hi = np.minimum(ka, kb), np.maximum(ka, kb)
    out = defaultdict(lambda: [0., 0., 0.])
    for m, j in ((sub & oo, 0), (sub & ~oo, 1)):
        for a, b, l, c in zip(lo[m], hi[m], L[m], cv[m]):
            out[(int(a), int(b))][j] += float(l)
            out[(int(a), int(b))][2] += float(l*c)
    return dict(out)


def main(dirs):
    g = None
    tot = defaultdict(float)
    tot_cs = defaultdict(float)
    onset_io = defaultdict(float)
    nwin = 0
    ev_count = defaultdict(int)
    for sd in map(Path, dirs):
        st = json.loads((sd/'stats.json').read_text())
        if g is None:
            g = CubeSphere(int(st['n']))
        parent = {}
        for e in st.get('event_list', []):
            if 'new_plate' in e:
                parent[int(e['new_plate'])] = (int(e['plate']), e['kind'], float(e['t']))
                ev_count[e['kind']] += 1
        start_ids = set(range(12))

        def line(k):                                  # ancestors of k, nearest first, with kinds
            out = []
            while k in parent:
                pk, kind, t = parent[k]
                out.append((pk, kind))
                k = pk
            return out

        def klass(a, b):
            if a in start_ids and b in start_ids:
                return 'start-start'
            for x, y in ((a, b), (b, a)):
                for anc, kind in line(x):
                    if anc == y:
                        return 'sibling:'+line(x)[0][1]
            young = max(a, b)
            return 'created:'+(parent[young][1] if young in parent else '?')
        snaps = {}
        for t in ALL_T:
            f = to.snap_path(sd, t)
            if Path(str(f)+'.npz').exists():
                snaps[t] = pair_lengths(g, f)
        times = sorted(snaps)
        for t in WIN_T:
            if t not in snaps:
                continue
            nwin += 1
            for k, (io, cs, lc) in snaps[t].items():
                c = klass(*k)
                tot[c] += io; tot_cs[c] += cs
                # onset: walk back while the pair has convergent contour
                i = times.index(t)
                while i > 0 and k in snaps[times[i-1]]:
                    i -= 1
                age = t-times[i]
                b = '<10' if age < 10 else '10-40' if age < 40 else '40-100' if age < 100 else '>=100'
                onset_io[b] += io
        print(sd.name, 'done', flush=True)
    s_io, s_cs = sum(tot.values()), sum(tot_cs.values())
    print(f"window -100 ... 0 (every 10 Myr), {len(dirs)} seeds: mean io {s_io/nwin/1e3:.1f}k km, cont-side {s_cs/nwin/1e3:.1f}k km")
    for c in sorted(tot, key=lambda c: -tot[c]):
        print(f"  {c:32s} io {tot[c]/nwin/1e3:6.1f}k ({tot[c]/s_io:.2f})   cont-side {tot_cs[c]/nwin/1e3:6.1f}k ({tot_cs[c]/max(s_cs,1):.2f})")
    print('  io by time since the pair began converging:', ', '.join(f"{b} Myr {onset_io[b]/nwin/1e3:.1f}k ({onset_io[b]/s_io:.2f})" for b in ('<10', '10-40', '40-100', '>=100')))
    print('  events per seed:', {k: v/len(dirs) for k, v in ev_count.items()})


if __name__ == '__main__':
    main(sys.argv[1:])
