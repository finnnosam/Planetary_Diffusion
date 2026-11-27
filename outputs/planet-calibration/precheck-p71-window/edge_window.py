"""Pass 71: K1 (passive edge), the edge's nearest-boundary split and the unbroken-rift edge (ridge-near, no own floor;
ridge_near2's class) read as window means over the saved snapshots -40 ... 0 Myr (every 5 Myr), not one 0 Myr snapshot.
Same instrument as edge_split.py / ridge_near2.py (CRUST1.0 1-degree cells, label-free, boundary within 300 km).
python edge_window.py <arm dirs...>"""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import to_crust1_cells
from planet_tectonics.graphops import nearest_label
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean())); k = lambda m: m.sum()*sp/1e6
T = ['t+0000.0' if t == 0 else f't-{-t:04d}.0' for t in range(-40, 1, 5)]
KEYS = ('edge', 'passive', 'ridge', 'unbroken', 'subduction', 'collision', 'transform')
def read(z, p):
    col = column_thickness_m(z, p); fl = z['continental'].astype(bool)
    cont = label_free_continent(to_crust1_cells(g, col, fl)[0])
    own = z['plate_id'].astype(np.int64); bt = z['boundary_type']
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    bl, bd = nearest_label(np.where(bt > 0, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    J = np.maximum(bl, 0); near = edge & (bd <= 300e3); ty = bt[J]
    el, ed = nearest_label(np.where(edge, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    oc = ~fl & (ed <= 300e3) & (el >= 0) & (own == own[np.maximum(el, 0)])
    n_oc = np.bincount(el[oc], minlength=g.size)
    ri = near & (ty == 1)
    return dict(edge=k(edge), passive=k(edge & ~near), ridge=k(ri), unbroken=k(ri & (n_oc == 0)), subduction=k(near & (ty == 2)),
                collision=k(near & ((ty == 3) | (ty == 5))), transform=k(near & (ty == 4)))
for arm in map(Path, sys.argv[1:]):
    S = []
    for s in sorted(arm.glob('seed*')):
        p = run_params(s); R = [read(np.load(s/'snapshots'/(t+'.npz')), p) for t in T]
        m = {q: np.mean([r[q] for r in R]) for q in KEYS}; S.append(m)
        print(f"{arm.name}/{s.name}: " + ' | '.join(f"{q} {m[q]:.1f}k" for q in KEYS) +
              " | unbroken by snapshot " + ' '.join('%.0f' % r['unbroken'] for r in R), flush=True)
    print(f"{arm.name} MEAN (window -40..0, 6 seeds): " + ' | '.join(f"{q} {np.mean([m[q] for m in S]):.1f}k" for q in KEYS) +
          f" | seed sd of unbroken {np.std([m['unbroken'] for m in S]):.1f}k, of passive {np.std([m['passive'] for m in S]):.1f}k", flush=True)
