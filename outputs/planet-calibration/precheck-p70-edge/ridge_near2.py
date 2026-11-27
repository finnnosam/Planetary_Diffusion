"""Pass 70: ridge_near.py's A = 0 class split: the edge's 300 km region holds no ocean floor of its own plate (the
ridge-typed boundary lies in continental crust: a rift before breakup, or a boundary along the coast), or young floor.
Also: the crust (km, crust + sediment) and continental flag at the nearest ridge node.  python ridge_near2.py <arm>"""
import sys, json, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import to_crust1_cells
from planet_tectonics.graphops import nearest_label
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean())); k = lambda m: m.sum()*sp/1e6
for s in sorted(Path(sys.argv[1]).glob('seed*')):
    z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s)
    col = column_thickness_m(z, p); fl = z['continental'].astype(bool)
    cont = label_free_continent(to_crust1_cells(g, col, fl)[0])
    own = z['plate_id'].astype(np.int64); bt = z['boundary_type']
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    bl, bd = nearest_label(np.where(bt > 0, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    J = np.maximum(bl, 0); ri = edge & (bd <= 300e3) & (bt[J] == 1)
    el, ed = nearest_label(np.where(edge, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    oc = ~fl & (ed <= 300e3) & (el >= 0) & (own == own[np.maximum(el, 0)])
    n_oc = np.bincount(el[oc], minlength=g.size)
    A = np.zeros(g.size); np.maximum.at(A, el[oc], z['ocean_age_myr'][oc].astype(float))
    none = ri & (n_oc == 0); young = ri & (n_oc > 0) & (A < 10)
    rc = col[J]/1e3; rf = fl[J]
    print(f"{s.name} ridge_near {k(ri):.1f}k | no own floor {k(none):.1f}k (ridge node flagged continental {k(none & rf):.1f}k, "
          f"its column < 20 km {k(none & (rc < 20)):.1f}k) | own floor < 10 Myr {k(young):.1f}k | d < 100 km: none {k(none & (bd < 1e5)):.1f}k "
          f"young {k(young & (bd < 1e5)):.1f}k | ridge node column p50: none {np.median(rc[none]) if none.any() else np.nan:.1f} km, young "
          f"{np.median(rc[young]) if young.any() else np.nan:.1f} km", flush=True)
