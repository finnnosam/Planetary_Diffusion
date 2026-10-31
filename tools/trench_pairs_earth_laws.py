"""Forty-fourth pass: do the plate layer's force laws, given Earth's plates and boundaries, reproduce Earth's
co-moving neighbours? rho < 0.5 share by Bird class and C(d), for Bird's observed poles, the torque law with PARAMS
and the length law (both from planet_reviews.force_balance_earth), all NNR. Run from the repo root with PYTHONPATH set."""
import numpy as np
from planet_reviews import trench_pairs as tp, force_balance_earth as fbe

r = fbe.evaluate_params(217)
g, pid, names, cont, st, K = r['g'], r['pid'], r['names'], r['cont'], r['st'], r['K']
E = tp.earth_segments()
oo = E['io'] | E['ocb']
x = st['x']; L = E['L']; cls = E['cls']; ms = fbe.R*1e3*fbe.YR
laws = {'observed (Bird NNR)': r['w_obs'], 'torque law, PARAMS': r['w'],
        'length law, PARAMS': fbe.remove_net_rotation(g, pid, fbe.length_law_omega(g, pid, cont, st, K))}
for lab, w in laws.items():
    va = np.cross(np.stack([w[int(k)] for k in st['A']]), x)*ms; vb = np.cross(np.stack([w[int(k)] for k in st['B']]), x)*ms
    nA = st['nA']
    close = -np.einsum('ij,ij->i', vb-va, nA)          # B relative to A along the normal out of A: + opening -> closing is minus
    rho = tp.rho(va, vb)
    row = [f"oo conv {tp.share(rho[oo] < .5, L[oo]):.2f} (closing p50 {tp.wq(close[oo], L[oo], (.5,))[0]:.0f}, <20: {tp.share(close[oo] < 20, L[oo]):.2f})"]
    for c in ('OSR', 'CRB', 'OTF', 'CTF', 'SUB', 'OCB', 'CCB'):
        m = cls == c; row.append(f'{c} {tp.share(rho[m] < .5, L[m]):.2f}')
    wn = np.stack([w[int(k)] for k in pid])
    coh = tp.coherence(g, pid, np.cross(wn, g.xyz), np.random.default_rng(0))
    print(f'{lab}\n  rho < 0.5 share: ' + '; '.join(row) + '\n  C(d): ' + tp.coh_txt(coh))
