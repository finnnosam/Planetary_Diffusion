import glob, sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
g = CubeSphere(217); A = g.area/g.area.sum()
for arm in sys.argv[1:]:
    fs = sorted(glob.glob(f'{arm}/seed1007/snapshots/t*.npz'), key=lambda s: float(s.split('/')[-1][1:-4]))
    worst = 0.
    for f in fs[-9:]:
        d = np.load(f); own = d['plate_id'] if 'plate_id' in d.files else d['owner']
        det = 0.
        for pid in np.unique(own):
            m = own == pid; lab = P.components(m, g.nbr); lab_m = lab[m]
            if lab_m.max() > 0:
                s = np.bincount(lab_m, weights=A[m]); det += s.sum()-s.max()
        worst = max(worst, det)
    print(arm, f'detached area max over the last 40 Myr {worst:.4f}')
