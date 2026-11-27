"""Pass 81 (journal 105): where the planet's L-T misfit sits (one arm, six seeds, a snapshot; default 0 Myr): per crust
bin (crust + sediment), area, the share of nodes on lithosphere at least 250 Myr old (at 0 Myr: untouched since the
start), and the area-weighted L mean of all / old / young nodes in the run's frame, beside Earth's Hoggard et al. 2020
SL2013sv node means. Run from the repo root: where_misfit.py <arm> [snapshot, e.g. t-0250.0]."""
import sys, numpy as np
from planet_tectonics.lithosphere import thickness_km
from planet_reviews.crust_review import run_params
from planet_grid.sphere_graph import CubeSphere
arm = sys.argv[1]; snap = sys.argv[2] if len(sys.argv) > 2 else 't+0000.0'
EARTH = {(20, 25): 97, (25, 30): 98, (30, 35): 107, (35, 40): 143, (40, 45): 158, (45, 55): 157, (55, 80): 197}
A = CubeSphere(217).area
TT, LL, AG, W = [], [], [], []
for s in range(1001, 1007):
    p = run_params(f'{arm}/seed{s}')
    q = np.load(f'{arm}/seed{s}/snapshots/{snap}.npz')
    la = q['lithosphere_age_myr'].astype(float); t = (q['crust_thickness_m']+q['sediment_m'])/1e3
    m = q['continental'].astype(bool) & np.isfinite(la) & (t >= 20)
    TT.append(t[m]); LL.append(thickness_km(la[m], p)); AG.append(la[m]); W.append(A[m])
TT, LL, AG, W = map(np.concatenate, (TT, LL, AG, W))
wm = lambda k: (LL[k]*W[k]).sum()/max(W[k].sum(), 1e-30)
print(f'{arm} {snap}: continental nodes {TT.size}; lithosphere >= 250 Myr old: area share {W[AG >= 250].sum()/W.sum():.2f}; '
      f'L mean {wm(np.ones(TT.size, bool)):.1f} km')
print('  T bin | area | old share | L mean all / old / young | Earth')
for (a, b), e in EARTH.items():
    m = (TT >= a) & (TT < b); o = m & (AG >= 250); y = m & (AG < 250)
    print(f'  {a}-{b} | {W[m].sum()/W.sum():.3f} | {W[o].sum()/max(W[m].sum(), 1e-30):.2f} | {wm(m):.0f} / {wm(o):.0f} / {wm(y):.0f} | {e}')
