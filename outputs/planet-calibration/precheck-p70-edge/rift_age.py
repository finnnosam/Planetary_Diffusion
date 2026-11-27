"""Pass 70: the planet's unbroken rifts (ridge-near edge with no own-plate floor within 300 km; ridge_near2.py): the
opening rate at the ridge node (full, mm/yr), the rift's age (oldest lithosphere age, C's clock, of the stretched crust
-- flagged continental, column < 20 km -- in the edge's 300 km region; C resets it when a gap is filled), the stretched
zone's width per edge node, and the crust (no sediment) at the ridge node.   python rift_age.py <arm>"""
import sys, numpy as np
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
    reg = (ed <= 300e3) & (el >= 0) & (own == own[np.maximum(el, 0)])
    n_oc = np.bincount(el[reg & ~fl], minlength=g.size)
    st = reg & fl & (col < 20e3)
    age = np.zeros(g.size); np.maximum.at(age, el[st], z['lithosphere_age_myr'][st].astype(float))
    wid = np.bincount(el[st], g.area[st], minlength=g.size)/sp/1e3          # km of stretched crust per edge node
    none = ri & (n_oc == 0)
    op = -z['convergence_mm_yr'][J][none]; a = age[none]; w = wid[none]; tc = z['crust_thickness_m'][J][none]/1e3
    q = lambda v: '/'.join(f'{x:.0f}' for x in np.percentile(v, [25, 50, 75]))
    print(f"{s.name} unbroken {k(none):.1f}k | opening p25/50/75 {q(op)} mm/yr | rift age {q(a)} Myr | stretched width "
          f"{q(w)} km | ridge-node crust {q(tc)} km | age > 25 Myr {k(none & (age > 25)):.1f}k, opening < 10 mm/yr "
          f"{k(none & (-z['convergence_mm_yr'][J] < 10)):.1f}k", flush=True)
