"""Pass 70: the ridge-near continental edge (non-passive, nearest boundary a ridge within 300 km; K1's instrument as
edge_cause.py) by how long its ocean has been opening: A = the oldest ocean floor of the edge's own plate within its
300 km neighbourhood (the floor between the edge and its ridge; not the downgoing floor, which lies on another plate),
the implied half-rate d / A, and whether the ridge's two plates are a rift pair (one split from the other, or both
from one parent) with the split's age.   python ridge_near.py <arm dirs...>"""
import sys, json, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import to_crust1_cells
from planet_tectonics.graphops import nearest_label

g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean()))
k = lambda m: m.sum()*sp/1e6


def one(s):
    z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s)
    ev = json.load(open(s/'stats.json'))['event_list']
    T1 = to_crust1_cells(g, column_thickness_m(z, p), z['continental'].astype(bool))[0]
    cont = label_free_continent(T1)
    own = z['plate_id'].astype(np.int64); bt = z['boundary_type']
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    bl, bd = nearest_label(np.where(bt > 0, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    near = edge & (bd <= 300e3); J = np.maximum(bl, 0)
    ri = near & (bt[J] == 1)
    # A: each edge node's Voronoi region within 300 km (graph distance), own plate, ocean floor (not continental)
    el, ed = nearest_label(np.where(edge, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    oc = ~z['continental'].astype(bool) & (ed <= 300e3) & (el >= 0) & (own == own[np.maximum(el, 0)])
    A = np.zeros(g.size); np.maximum.at(A, el[oc], z['ocean_age_myr'][oc].astype(float))
    # rift pairs: parent -> child edges of rift events; sisters share a parent too
    par = {}; tsplit = {}
    for e in ev:
        if e['kind'] == 'rift':
            par[e['new_plate']] = e['plate']; tsplit[e['new_plate']] = e['t']
    def pair_age(a, b):
        if par.get(a) == b: return -tsplit[a]
        if par.get(b) == a: return -tsplit[b]
        if a in par and par.get(a) == par.get(b): return -max(tsplit[a], tsplit[b])
        return np.nan
    I = np.flatnonzero(ri)
    oth = np.full(g.size, -1, np.int64)
    for i in np.flatnonzero(bt == 1):
        js = [j for j in nb[i] if j >= 0 and own[j] != own[i]]
        if js: oth[i] = max(set(own[js]), key=list(own[js]).count)
    pa = np.array([pair_age(int(own[J[i]]), int(oth[J[i]])) for i in I])
    a_, d_ = A[I], bd[I]/1e3
    r = {'ridge_near': k(ri)}
    for lo, hi in ((0, 10), (10, 25), (25, 50), (50, 100), (100, 1e9)):
        m = (a_ >= lo) & (a_ < hi)
        r[f'A{lo}-{int(min(hi, 999))}'] = m.sum()*sp/1e6
        r[f'hr{lo}'] = float(np.median(d_[m]/np.maximum(a_[m], .5)*.1)) if m.any() else np.nan   # cm/yr
    r['pair_rift'] = np.isfinite(pa).sum()*sp/1e6
    r['pair_rift_lt25'] = (pa < 25).sum()*sp/1e6
    r['pair_none'] = (~np.isfinite(pa)).sum()*sp/1e6
    r['pair_none_A_ge50'] = ((~np.isfinite(pa)) & (a_ >= 50)).sum()*sp/1e6
    return r


out = {}
for arm in map(Path, sys.argv[1:]):
    rs = []
    for s in sorted(arm.glob('seed*')):
        r = one(s); rs.append(r)
        print(f'  {s.name} ' + ' '.join(f'{x} {v:.1f}' for x, v in r.items()), flush=True)
    m = {x: float(np.nanmean([r[x] for r in rs])) for x in rs[0]}
    print(f'{arm.name} mean ' + ' '.join(f'{x} {v:.1f}' for x, v in m.items()), flush=True)
    out[arm.name] = {'seeds': rs, 'mean': m}
json.dump(out, open(sys.argv[1].rstrip('/').split('/')[-1]+'-ridge_near.json', 'w'), indent=1, default=float)
