"""Bit-identity fingerprint for speed work (journal 154): n, seed, threads T, K steps -> one md5 per step of every
state() array, every plate's rasters and rotation, the mantle flow record and the ledgers (erosion_totals).
  python outputs/planet-calibration/precheck-speed/ident.py 48 7 2 70 out.json
Two fingerprints from the same arguments must be equal line by line (compare.py)."""
import sys, hashlib, json, time, glob
sys.path[:0] = glob.glob('packages/*/src')
import numpy as np
n, seed, T, K, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim


def md5(a):
    return hashlib.md5(np.ascontiguousarray(a).tobytes()).hexdigest()


g = CubeSphere(n)
sim = PlateSim(g, seed=seed, t0=-250., params={'threads': T})
rows = []


def rec(s):
    st = s.state()
    r = {k: md5(v) for k, v in st.items() if isinstance(v, np.ndarray)}
    for pl in s.plates:
        for k, v in vars(pl).items():
            if isinstance(v, np.ndarray):
                r[f'plate{pl.id}.{k}'] = md5(v)
    if getattr(s, 'mantle', None) is not None and s.mantle.flow is not None:
        r['_flowH'] = md5(s.mantle.flow.H)
    r['_totals'] = hashlib.md5(json.dumps(s.erosion_totals, sort_keys=True, default=float).encode()).hexdigest()
    r['_events'] = len(s.events)
    rows.append(r)


t0 = time.perf_counter()
rec(sim)
sim.run(K, rec, 1)
json.dump({'seconds': time.perf_counter()-t0, 'rows': rows}, open(out, 'w'))
print('done', len(rows), f'{time.perf_counter()-t0:.1f} s')
