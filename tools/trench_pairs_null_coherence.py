"""Forty-fourth pass: C(d) of different plates under random pole directions (each plate keeps its speed), Earth and
planet seeds; the null for planet_reviews.trench_pairs. Args: seed folders. Run from the repo root with PYTHONPATH set."""
import sys, json, numpy as np
from planet_reviews import trench_pairs as tp, force_balance_earth as fbe
from planet_grid.sphere_graph import CubeSphere
def null(g, pid, w_obs, n=12, seed=0):
    rng = np.random.default_rng(seed); out = []
    ids = np.unique(pid)
    mags = [np.linalg.norm(w_obs[int(k)]) for k in ids]
    for r in range(n):
        w = {int(k): rng.normal(size=3) for k in ids}
        w = {k: v/np.linalg.norm(v)*rng.choice(mags) for k, v in w.items()}
        w = fbe.remove_net_rotation(g, pid, w)
        wn = np.stack([w[int(k)] for k in pid])
        out.append(tp.coherence(g, pid, np.cross(wn, g.xyz), rng, n_ref=400, n_cmp=20000))
    return np.mean(out, 0), np.std(out, 0)
def neff(g, pid):
    a = np.array([g.area[pid == k].sum() for k in np.unique(pid)]); return a.sum()**2/(a**2).sum()
g, pid, names, _ = fbe.earth_plates(217)
w = fbe.earth_omega(g, pid, names)
m, s = null(g, pid, w)
print('Earth N_eff %.1f; random poles C(d): %s (sd %s)' % (neff(g, pid), tp.coh_txt(m), ' '.join('%.2f' % x for x in s)))
for sd in sys.argv[1:]:
    z = np.load(sd+'/snapshots/t+0000.0.npz'); pl = json.load(open(sd+'/snapshots/t+0000.0-plates.json'))
    pid = z['plate_id'].astype(np.int64)
    w = {int(q['id']): np.asarray(q['axis'])*q['speed_cm_yr'] for q in pl}
    m, s = null(g, pid, w, n=6)
    print('%s N_eff %.1f; random poles C(d): %s' % (sd[-8:], neff(g, pid), tp.coh_txt(m)))
