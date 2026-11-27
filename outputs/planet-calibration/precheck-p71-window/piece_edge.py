"""Pass 71: where the continental edge lies by the size of its continental piece (connected label-free continent on
n = 217, CRUST1.0 1-degree instrument), passive or not; Earth (CRUST1.0, PB2002) against planet arms (window -40..0 Myr).
python piece_edge.py <arm dirs...>"""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_grid.paths import DATA
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import earth_nodes, to_crust1_cells
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_tectonics.graphops import nearest_label
g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))
BINS = [0, 0.1, 1, 10, 1e9]; NAMES = ['<0.1', '0.1-1', '1-10', '>10']
def comps(m):
    lab = np.where(m, np.arange(g.size), -1)
    while True:
        nl = lab.copy()
        for j in range(nb.shape[1]):
            q = np.where(ok[:, j], nb[:, j], 0); v = np.where(ok[:, j] & m & m[q], lab[q], -1); nl = np.maximum(nl, v)
        nl = np.where(m, nl, -1)
        if (nl == lab).all(): return lab
        lab = nl
def read(cont, isb):
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    _, bd = nearest_label(np.where(isb, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    pas = edge & (bd > 300e3)
    lab = comps(cont); A = np.bincount(np.maximum(lab, 0), weights=np.where(cont, a, 0), minlength=g.size)/1e12
    sz = A[np.maximum(lab, 0)]; r = {}
    for i, n in enumerate(NAMES):
        s = cont & (sz >= BINS[i]) & (sz < BINS[i+1])
        r[n] = (s.sum() and a[s].sum()/1e12, (edge & s).sum()*sp/1e6, (pas & s).sum()*sp/1e6, len(np.unique(lab[s])))
    return r
def fmt(nm, r):
    return nm + ': ' + ' | '.join(f"{n} Mkm2: {r[n][3]:.0f} pieces {r[n][0]:.1f} Mkm2 edge {r[n][1]:.1f}k passive {r[n][2]:.1f}k" for n in NAMES)
T, _ = earth_nodes(g); s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
x = xyz_from_lonlat(s['STARTLONG'], s['STARTLAT'])+xyz_from_lonlat(s['FINALLONG'], s['FINALLAT']); x /= np.linalg.norm(x, axis=1, keepdims=True)
isb = np.zeros(g.size, bool); isb[fields.nearest_nodes(g, x)] = True
if __name__ == '__main__':
    print(fmt('Earth', read(label_free_continent(T), isb)), flush=True)
    for arm in map(Path, sys.argv[1:]):
        R = []
        for sd in sorted(arm.glob('seed*')):
            p = run_params(sd)
            for t in ['t+0000.0'] + [f't-{q:04d}.0' for q in (10, 20, 30, 40)]:
                z = np.load(sd/'snapshots'/(t+'.npz'))
                R.append(read(label_free_continent(to_crust1_cells(g, column_thickness_m(z, p), z['continental'].astype(bool))[0]), z['boundary_type'] > 0))
        m = {n: tuple(np.mean([r[n][i] for r in R]) for i in range(4)) for n in NAMES}
        print(fmt(f'{arm.name} (6 seeds x 5 snapshots, -40..0)', m), flush=True)
