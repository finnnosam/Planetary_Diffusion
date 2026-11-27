# Pass 67: for planet continental edge with margin age >= 50 Myr that is NOT passive (a boundary within 300 km):
# which boundary type is nearest (share of length). Also the plates' speeds: continental-edge plates vs the rest.
import sys, json, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_tectonics.graphops import nearest_label
from planet_tectonics.params import RIDGE, SUBDUCTION, COLLISION, TRANSFORM
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean()))
NAMES = {RIDGE: 'ridge', SUBDUCTION: 'subduction', COLLISION: 'collision', TRANSFORM: 'transform'}
acc = {}; accy = {}
for s in map(Path, sys.argv[1:]):
    z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s)
    cont = label_free_continent(column_thickness_m(z, p)); bt = z['boundary_type']
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    age = np.where(cont, np.nan, z['ocean_age_myr'].astype(float))
    lab, d = nearest_label(np.where(edge, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    sel = (~cont) & np.isfinite(age) & (d <= 300e3) & (lab >= 0)
    m = np.full(g.size, np.nan); np.fmax.at(m, lab[sel], age[sel])
    bl, bd = nearest_label(np.where(bt != 0, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    for tgt, lo, hi in ((acc, 50, 1e9), (accy, 0, 50)):
        f = edge & (m >= lo) & (m < hi) & (bd <= 300e3)
        for k, v in NAMES.items(): tgt[v] = tgt.get(v, 0) + (bt[bl[f]] == k).sum()*sp/1e3
        tgt['n'] = tgt.get('n', 0) + 1
for nm, t in (('margin age >= 50 Myr', acc), ('margin age < 50 Myr', accy)):
    tot = sum(v for k, v in t.items() if k != 'n')
    print(f'{Path(sys.argv[1]).parent.name} not passive, {nm}: {tot/t["n"]/1e3:.1f}k km per seed; nearest boundary ' +
          ', '.join(f'{k} {v/tot:.2f}' for k, v in t.items() if k != 'n'))
