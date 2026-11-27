"""Pass 69: K1's two factors on the planet and on Earth, same instrument (CRUST1.0 1-degree cells, label-free, nodes of
n = 217): the label-free continental edge per Mkm^2 of continent, its passive share (no plate boundary within 300 km),
and the non-passive edge by its nearest boundary's type.   python edge_split.py <arm dirs...>"""
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
TYPES = ('ridge', 'subduction', 'collision', 'transform')
PLANET = {1: 'ridge', 2: 'subduction', 3: 'collision', 4: 'transform', 5: 'collision'}
EARTH = {'OSR': 'ridge', 'CRB': 'ridge', 'SUB': 'subduction', 'OCB': 'collision', 'CCB': 'collision', 'OTF': 'transform',
         'CTF': 'transform'}


def split(cont, btype):
    """btype: node array of type names ('' off boundary)."""
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    bl, bd = nearest_label(np.where(btype != '', np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    near = edge & (bd <= 300e3)
    L = edge.sum()*sp/1e6; A = a[cont].sum()/1e12
    r = {'edge_k_km': L, 'cont_mkm2': A, 'edge_per_mkm2': L*1e3/A, 'passive_k_km': (edge & ~near).sum()*sp/1e6,
         'passive_share': (edge & ~near).sum()/edge.sum()}
    for t in TYPES:
        r['near_'+t] = (near & (btype[np.maximum(bl, 0)] == t)).sum()*sp/1e6
    # nearest-boundary distance distribution of the non-passive edge (what a wider ocean would free)
    d = bd[near]/1e3
    r['near_lt100'] = (d < 100).sum()*sp/1e6; r['near_100_300'] = (d >= 100).sum()*sp/1e6
    return r


def fmt(nm, r):
    return (f"{nm:26s} edge {r['edge_k_km']:6.1f}k km on {r['cont_mkm2']:5.1f} Mkm2 ({r['edge_per_mkm2']:4.0f} km/Mkm2) | passive "
            f"{r['passive_k_km']:6.1f}k ({r['passive_share']:.2f}) | not passive, nearest " +
            ' / '.join(f"{t} {r['near_'+t]:5.1f}k" for t in TYPES) +
            f" | boundary < 100 km {r['near_lt100']:5.1f}k, 100-300 km {r['near_100_300']:5.1f}k")


def earth():
    T, _ = earth_nodes(g)
    s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
    x = xyz_from_lonlat(s['STARTLONG'], s['STARTLAT'])+xyz_from_lonlat(s['FINALLONG'], s['FINALLAT'])
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    bt = np.full(g.size, '', dtype=object)
    bt[fields.nearest_nodes(g, x)] = [EARTH.get(c.strip(), 'transform') for c in s['STEPCLASS']]
    return split(label_free_continent(T), bt.astype(str))


R = {'Earth': [earth()]}
print(fmt('Earth (CRUST1.0, PB2002)', R['Earth'][0]), flush=True)
for arm in map(Path, sys.argv[1:]):
    R[arm.name] = []
    for s in sorted(arm.glob('seed*')):
        z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s)
        T1 = to_crust1_cells(g, column_thickness_m(z, p), z['continental'].astype(bool))[0]
        bt = np.array([PLANET.get(int(v), '') for v in z['boundary_type']])
        r = split(label_free_continent(T1), bt); R[arm.name].append(r)
        print(fmt(f'  {arm.name[12:]}/{s.name}', r), flush=True)
    m = {k: float(np.mean([r[k] for r in R[arm.name]])) for k in R[arm.name][0]}
    print(fmt(f'{arm.name} mean', m), flush=True)
