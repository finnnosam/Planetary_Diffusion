"""Pass 81: the planet's young lithosphere (age < 250 Myr) under thick crust (crust + sediment > 45 km) at 0 Myr and
-5 Myr: how young (reset within 1 / 5 / 20 Myr), how thick the crust (> 65 km: the foundering root), and the area
share of continental thick crust it makes. Run from the repo root."""
import sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
A = CubeSphere(217).area
arm = sys.argv[1]
for snap in ('t+0000.0', 't-0005.0'):
    rows = []
    for s in range(1001, 1007):
        q = np.load(f'{arm}/seed{s}/snapshots/{snap}.npz')
        t = (q['crust_thickness_m']+q['sediment_m'])/1e3; la = q['lithosphere_age_myr'].astype(float)
        c = q['continental'].astype(bool) & np.isfinite(la)
        thick = c & (t > 45); y = thick & (la < 250)
        rows.append([A[thick].sum(), A[y].sum(), A[y & (la < 1)].sum(), A[y & (la < 5)].sum(), A[y & (la < 20)].sum(),
                     A[y & (t > 65)].sum(), A[thick & (t > 65)].sum(), A[y & (q['tectonic_thickening_m'] > 0)].sum()])
    r = np.array(rows).sum(0)
    print(f'{snap}: thick crust (> 45 km) area {r[0]/6e12:.2f} Mkm2/seed; young share {r[1]/r[0]:.2f}; of the young: age < 1 Myr '
          f'{r[2]/r[1]:.2f}, < 5 {r[3]/r[1]:.2f}, < 20 {r[4]/r[1]:.2f}, crust > 65 km {r[5]/r[1]:.2f}, thickening now {r[7]/r[1]:.2f}; '
          f'> 65 km share of all thick {r[6]/r[0]:.3f}')
