"""Probe (pass 42): cumulative convergence of the subduction zone at each place along the boundary (not per plate
pair). For 1 Myr snapshots -40 ... 0 Myr, each convergent contour segment (typed subduction) takes the cumulative
convergence of the nearest convergent segment of the previous snapshot within R rings of grid nodes (the boundary
moves with the overriding plate, <= ~80 km/Myr) and adds its own closing rate x 1 Myr; with none in reach it starts
from its own rate (a new subduction zone). Segments convergent at -40 Myr start as 'old'. Prints, for the window
-10 ... 0, the intra-oceanic and continental-side trench whose local cumulative convergence is below D km.

  python local_slab_probe.py <rings> <seed dir> [...]
"""
import json
import sys
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import boundary_lengths as bl
from planet_reviews import trench_origin as to

D_KM = (50., 100., 125., 150., 250.)
TIMES = [float(t) for t in range(-40, 1)]
OLD = 1e9


def segs(g, f):
    z = np.load(str(f)+'.npz'); plates = json.loads(Path(str(f)+'-plates.json').read_text())
    pid = z['plate_id']; cont = z['continental'].astype(bool)
    s = bl.typed_segments(g, pid, plates, cont)
    sub = s['typ'] == 2
    oo = ~(cont[s['na']] | cont[s['nb']])
    return s['L'][sub], s['conv'][sub], oo[sub], fields.nearest_nodes(g, s['mid'][sub])


def dilate(g, F, rings):
    nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
    for _ in range(rings):
        F = np.maximum(F, F[nb].max(1))
    return F


def main(rings, dirs):
    g = None
    out = []
    for sd in map(Path, dirs):
        st = json.loads((sd/'stats.json').read_text())
        if g is None:
            g = CubeSphere(int(st['n']))
        F = None
        rows = []
        for i, t in enumerate(TIMES):
            L, cv, oo, nn = segs(g, to.snap_path(sd, t))
            if F is None:
                cum = np.full(L.size, OLD)
            else:
                prev = F[nn]
                cum = np.where(prev >= 0, prev, 0.)+cv          # mm/yr x 1 Myr = km
                cum = np.where(prev >= OLD, OLD, cum)
            G = np.full(g.size, -1.)
            np.maximum.at(G, nn, cum)
            F = dilate(g, G, rings)
            if t >= -10:
                r = {'io': L[oo].sum(), 'cs': L[~oo].sum(), 'io_old': L[oo & (cum >= OLD)].sum()}
                for D in D_KM:
                    r[f'io{D:.0f}'] = L[oo & (cum < D)].sum(); r[f'cs{D:.0f}'] = L[~oo & (cum < D)].sum()
                rows.append(r)
        m = {k: np.mean([r[k] for r in rows]) for k in rows[0]}
        out.append(m)
        print(f"{sd.name}: io {m['io']/1e3:.1f}k (old >= 40 Myr {m['io_old']/1e3:.1f}k), cs {m['cs']/1e3:.1f}k; "
              + '; '.join(f"< {D:.0f}: io {m[f'io{D:.0f}']/1e3:.1f}k cs {m[f'cs{D:.0f}']/1e3:.1f}k" for D in D_KM), flush=True)
    M = {k: np.mean([r[k] for r in out]) for k in out[0]}
    print(f"mean ({rings} rings): io {M['io']/1e3:.1f}k (old {M['io_old']/1e3:.1f}k), cs {M['cs']/1e3:.1f}k; "
          + '; '.join(f"< {D:.0f}: io {M[f'io{D:.0f}']/1e3:.1f}k cs {M[f'cs{D:.0f}']/1e3:.1f}k" for D in D_KM))


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2:])
