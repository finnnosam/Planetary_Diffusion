"""Pass 70: the unbroken rifts (ridge_near2.py's class) by their plate pair: time since the split that made the pair
(event list; a pair is parent/child or two children of one parent), the edge km per pair, the pair's opening (full rate
at the ridge nodes, mm/yr) and the column (crust + sediment) along the pair's ridge.   python rift_pairs.py <arm>"""
import sys, json, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import to_crust1_cells
from planet_tectonics.graphops import nearest_label
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean()))
for s in sorted(Path(sys.argv[1]).glob('seed*')):
    z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s); ev = json.load(open(s/'stats.json'))['event_list']
    par = {}; ts = {}
    for e in ev:
        if e['kind'] == 'rift': par[e['new_plate']] = e['plate']; ts[e['new_plate']] = e['t']
    def age(a, b):
        if par.get(a) == b: return -ts[a]
        if par.get(b) == a: return -ts[b]
        if a in par and par.get(a) == par.get(b): return -max(ts[a], ts[b])
        return -1.
    col = column_thickness_m(z, p); fl = z['continental'].astype(bool)
    cont = label_free_continent(to_crust1_cells(g, col, fl)[0])
    own = z['plate_id'].astype(np.int64); bt = z['boundary_type']
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    bl, bd = nearest_label(np.where(bt > 0, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    J = np.maximum(bl, 0); ri = edge & (bd <= 300e3) & (bt[J] == 1)
    el, ed = nearest_label(np.where(edge, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    reg = (ed <= 300e3) & (el >= 0) & (own == own[np.maximum(el, 0)])
    none = ri & (np.bincount(el[reg & ~fl], minlength=g.size) == 0)
    oth = np.full(g.size, -1)
    for i in np.flatnonzero(bt == 1):
        js = [j for j in nb[i] if j >= 0 and own[j] != own[i]]
        if js: oth[i] = max(set(own[js]), key=list(own[js]).count)
    pairs = {}
    for i in np.flatnonzero(none):
        key = tuple(sorted((int(own[J[i]]), int(oth[J[i]])))); pairs.setdefault(key, []).append(i)
    rows = []
    for key, I in pairs.items():
        rnodes = (bt == 1) & (((own == key[0]) & (oth == key[1])) | ((own == key[1]) & (oth == key[0])))
        rows.append((len(I)*sp/1e3, age(*key), np.median(-z['convergence_mm_yr'][rnodes]) if rnodes.any() else np.nan,
                     np.median(col[rnodes])/1e3 if rnodes.any() else np.nan, (fl & rnodes).sum()/max(rnodes.sum(), 1), key))
    rows.sort(reverse=True)
    tot = sum(r[0] for r in rows)
    print(f"{s.name} unbroken {tot/1e3:.1f}k km in {len(rows)} pairs; by split age: none (not a rift pair) "
          f"{sum(r[0] for r in rows if r[1] < 0)/1e3:.1f}k, < 25 Myr {sum(r[0] for r in rows if 0 <= r[1] < 25)/1e3:.1f}k, "
          f"25-50 {sum(r[0] for r in rows if 25 <= r[1] < 50)/1e3:.1f}k, >= 50 {sum(r[0] for r in rows if r[1] >= 50)/1e3:.1f}k", flush=True)
    for r in rows[:5]:
        print(f"    pair {r[5]}: {r[0]:6.0f} km, split {r[1]:5.0f} Myr ago, opening p50 {r[2]:4.0f} mm/yr, ridge column p50 "
              f"{r[3]:4.1f} km, ridge nodes flagged continental {r[4]:.2f}", flush=True)
