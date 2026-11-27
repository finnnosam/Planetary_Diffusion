"""Pass 70: the choosing reading (journal 93). The non-passive continental edge (K1's instrument, as edge_split.py of
pass 69: CRUST1.0 1-degree cells, label-free, n = 217, boundary within 300 km) by the motion at its nearest boundary
node: the two plates' relative velocity there (normal = the snapshot's convergence rate, shear = the rest), and the
nearest boundary's type. Transform-typed nodes split into 'stalled' (|v_rel| < 5 mm/yr: the typing's vmin; a pair that
has stopped moving apart) and 'oblique' (real shear). Earth: PB2002 steps with their own velocities (VELOCITYDI,
VELOCITYRI, mm/yr), by Bird's class and retyped by the planet's rule (|normal| < max(5, shear/2) -> transform).
    python edge_cause.py <arm dirs...>"""
import sys, json, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_grid.paths import DATA
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import earth_nodes, to_crust1_cells
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_tectonics.graphops import nearest_label

g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))
PLANET = {1: 'ridge', 2: 'subduction', 3: 'collision', 4: 'transform', 5: 'collision'}
EARTH = {'OSR': 'ridge', 'CRB': 'ridge', 'SUB': 'subduction', 'OCB': 'collision', 'CCB': 'collision', 'OTF': 'transform',
         'CTF': 'transform'}
VMIN = 5.


def retype(vn, vt):
    """The planet's typing rule on a normal (+ closing) and shear rate, mm/yr; 'convergent' for closing."""
    return np.where(np.abs(vn) < np.maximum(VMIN, .5*np.abs(vt)), 'transform', np.where(vn > 0, 'convergent', 'ridge'))


def read(cont, btype, vn, vt):
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    bl, bd = nearest_label(np.where(btype != '', np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    near = edge & (bd <= 300e3)
    k = lambda m: m.sum()*sp/1e6
    J = np.maximum(bl, 0)
    t, n_, s_ = btype[J], vn[J], vt[J]
    v = np.hypot(n_, s_)
    r = {'edge': k(edge), 'passive': k(edge & ~near), 'near': k(near)}
    for ty in ('ridge', 'subduction', 'collision', 'transform'):
        r[ty] = k(near & (t == ty))
    tr = near & (t == 'transform')
    r['tr_stalled'] = k(tr & (v < VMIN)); r['tr_oblique'] = k(tr & (v >= VMIN))
    ri = near & (t == 'ridge')
    op = -n_[ri]
    r['ridge_open_p25'], r['ridge_open_p50'], r['ridge_open_p75'] = (np.percentile(op, [25, 50, 75]) if op.size else [np.nan]*3)
    r['ridge_lt10'] = k(ri & (-n_ < 10)); r['ridge_10_30'] = k(ri & (-n_ >= 10) & (-n_ < 30)); r['ridge_ge30'] = k(ri & (-n_ >= 30))
    r['ridge_d_lt100'] = k(ri & (bd < 100e3))
    # retyped by the planet's rule (Earth: from its own velocities; planet: identity check of the instrument)
    rt = retype(n_, s_)
    r['rule_transform'] = k(near & (rt == 'transform')); r['rule_ridge'] = k(near & (rt == 'ridge'))
    r['rule_tr_stalled'] = k(near & (rt == 'transform') & (v < VMIN))
    return r


def earth():
    T, _ = earth_nodes(g)
    s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
    x = xyz_from_lonlat(s['STARTLONG'], s['STARTLAT'])+xyz_from_lonlat(s['FINALLONG'], s['FINALLAT'])
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    i = fields.nearest_nodes(g, x)
    bt = np.full(g.size, '', dtype=object); vn = np.zeros(g.size); vt = np.zeros(g.size)
    bt[i] = [EARTH.get(c.strip(), 'transform') for c in s['STEPCLASS']]
    vn[i] = -np.asarray(s['VELOCITYDI'], float); vt[i] = np.asarray(s['VELOCITYRI'], float)
    return read(label_free_continent(T), bt.astype(str), vn, vt)


def planet(z, plates):
    om = {int(q['id']): np.array(q['axis'])*q['speed_cm_yr'] for q in plates}   # cm/yr at 90 degrees
    own = z['plate_id'].astype(np.int64); x = g.xyz
    bt = np.array([PLANET.get(int(v), '') for v in z['boundary_type']])
    vn = z['convergence_mm_yr'].astype(float).copy(); vt = np.zeros(g.size)
    for i in np.flatnonzero(bt != ''):
        js = [j for j in nb[i] if j >= 0 and own[j] != own[i]]
        if not js:
            continue
        b = max(set(own[js]), key=list(own[js]).count)
        dv = np.cross(om.get(own[i], np.zeros(3))-om.get(b, np.zeros(3)), x[i])*10.     # mm/yr
        vt[i] = np.sqrt(max(dv@dv-vn[i]**2, 0.))
    return bt, vn, vt


KEYS = ('edge', 'passive', 'near', 'ridge', 'subduction', 'collision', 'transform', 'tr_stalled', 'tr_oblique',
        'ridge_lt10', 'ridge_10_30', 'ridge_ge30', 'ridge_d_lt100', 'rule_transform', 'rule_ridge', 'rule_tr_stalled')


def fmt(nm, r):
    return (f"{nm:24s} " + ' '.join(f"{k} {r[k]:5.1f}k" for k in KEYS) +
            f" | ridge opening p25/50/75 {r['ridge_open_p25']:.0f}/{r['ridge_open_p50']:.0f}/{r['ridge_open_p75']:.0f} mm/yr")


R = {'Earth': earth()}
print(fmt('Earth', R['Earth']), flush=True)
out = {'Earth': R['Earth']}
for arm in map(Path, sys.argv[1:]):
    rs = []
    for s in sorted(arm.glob('seed*')):
        z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s)
        pl = json.load(open(s/'snapshots'/'t+0000.0-plates.json'))
        T1 = to_crust1_cells(g, column_thickness_m(z, p), z['continental'].astype(bool))[0]
        r = read(label_free_continent(T1), *planet(z, pl)); rs.append(r)
        print(fmt(f'  {s.name}', r), flush=True)
    m = {k: float(np.nanmean([r[k] for r in rs])) for k in rs[0]}
    print(fmt(f'{arm.name[:20]} mean', m), flush=True)
    out[arm.name] = {'seeds': rs, 'mean': m}
rv = Path(sys.argv[1]).resolve().parent/'reviews'          # beside the sweep's other reviews
rv.mkdir(exist_ok=True)
json.dump(out, open(rv/'edge_cause.json', 'w'), indent=1, default=float)
