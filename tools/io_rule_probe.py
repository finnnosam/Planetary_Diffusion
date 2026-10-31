"""Probe (pass 42): the planet's intra-oceanic trench measured by Earth's rule (trench_origin.earth_reference: a
subduction step is intra-oceanic when the crust is oceanic at d/2 and d to both sides of the step, d = 200 km; 100 km
beside) instead of the review's 0.75 cell. Planet crust is 'oceanic' where the continental flag is False. Window
-100 ... 0 Myr every 10 Myr, as T2.

  python io_rule_probe.py <seed dir> [...]
"""
import json
import sys
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import boundary_lengths as bl
from planet_reviews import trench_origin as to

WIN_T = [float(t) for t in range(-100, 1, 10)]


def rows(g, f):
    z = np.load(str(f)+'.npz'); plates = json.loads(Path(str(f)+'-plates.json').read_text())
    pid = z['plate_id']; cont = z['continental'].astype(bool)
    s = bl.typed_segments(g, pid, plates, cont)
    L, typ, mid, na, nb = s['L'], s['typ'], s['mid'], s['na'], s['nb']
    seg = bl.segments(g, pid)
    tang = seg[:, 1]-seg[:, 0]
    nrm = np.cross(mid, tang); nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-15)
    sub = typ == 2
    near = ~(cont[na] | cont[nb])
    out = {'near': float(L[sub & near].sum()), 'cs_near': float(L[sub & ~near].sum())}
    for d in (100., 200.):
        ok = np.ones(len(L), bool)
        for sg in (1., -1.):
            for x in (d/2, d):
                q = mid+sg*x*1e3/g.radius*nrm; q /= np.linalg.norm(q, axis=1, keepdims=True)
                ok &= ~cont[fields.nearest_nodes(g, q)]
        out[f'io{d:.0f}'] = float(L[sub & ok].sum())
        out[f'io{d:.0f}_and_near'] = float(L[sub & ok & near].sum())
    out['sub'] = float(L[sub].sum())
    return out


def main(dirs):
    g = None
    allr = []
    for sd in map(Path, dirs):
        st = json.loads((sd/'stats.json').read_text())
        if g is None:
            g = CubeSphere(int(st['n']))
        R = [rows(g, to.snap_path(sd, t)) for t in WIN_T if Path(str(to.snap_path(sd, t))+'.npz').exists()]
        m = {k: np.mean([r[k] for r in R]) for k in R[0]}
        allr.append(m)
        print(f"{sd.name}: subduction {m['sub']/1e3:.1f}k; io 0.75 cell (T2) {m['near']/1e3:.1f}k; Earth rule 200 km "
              f"{m['io200']/1e3:.1f}k, 100 km {m['io100']/1e3:.1f}k; continental side by Earth rule "
              f"{(m['sub']-m['io200'])/1e3:.1f}k", flush=True)
    M = {k: np.mean([r[k] for r in allr]) for k in allr[0]}
    print(f"mean: subduction {M['sub']/1e3:.1f}k; T2 rule {M['near']/1e3:.1f}k (cont side {M['cs_near']/1e3:.1f}k); "
          f"Earth rule io 200 km {M['io200']/1e3:.1f}k (cont side {(M['sub']-M['io200'])/1e3:.1f}k), 100 km "
          f"{M['io100']/1e3:.1f}k (cont side {(M['sub']-M['io100'])/1e3:.1f}k). Earth: 20.8k / 30.7k (200 km), 25.8k (100 km)")


if __name__ == '__main__':
    main(sys.argv[1:])
