"""Pass 95 pre-check reader: by total-thickness bin (crust + bulk sediment, flagged continent, as level_split bins the
planet), area means at 0 Myr and over the -40 ... 0 window: bulk sediment (m), elevation on the datum (m), share under
the datum, and the rock share of the sediment (sediment_rock 'burial'). Earth's means from level_split (pass 95 header)."""
import glob, sys
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.isostasy import sediment_column
from planet_reviews.crust_review import run_params

BINS = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55)]
EARTH_SED = [2607, 2777, 1876, 1317, 1616, 1076]
EARTH_Z = [-750, -196, 219, 441, 520, 1127]
g = CubeSphere(217); a = g.area
for arm in sys.argv[1:]:
    run = glob.glob(f'{arm}/seed*')[0]; p = run_params(run)
    fs = sorted(glob.glob(f'{run}/snapshots/t*.npz'), key=lambda s: float(s.split('/')[-1][1:-4]))
    fs = [f for f in fs if float(f.split('/')[-1][1:-4]) >= -40]
    acc = []
    for f in fs:
        d = np.load(f); sl = float(d['sea_level_m']); zd = d['elevation_m'].astype(float)+sl
        cont = d['continental'].astype(bool); sh = sediment_column(d['sediment_m'].astype(float), cont, p)[0]
        tot = (d['crust_thickness_m']+sh)/1e3
        rk = d['sediment_rock_m'].astype(float) if 'sediment_rock_m' in d.files else np.zeros(g.size)
        row = []
        for lo, hi in BINS:
            m = cont & (tot >= lo) & (tot < hi); w = a[m]
            row.append([(sh[m]*w).sum()/w.sum(), (zd[m]*w).sum()/w.sum(), w[zd[m] < 0].sum()/w.sum(),
                        (rk[m]*w).sum()/max((d['sediment_m'][m]*w).sum(), 1.), w.sum()/a.sum()])
        acc.append(row)
    acc = np.array(acc)
    print(arm.split('/')[-1])
    for nm, v in (('0 Myr', acc[-1]), ('window', acc.mean(0))):
        print(f'  {nm}:')
        for (lo, hi), r, es, ez in zip(BINS, v, EARTH_SED, EARTH_Z):
            print(f'    {lo}-{hi} km  share {r[4]:.3f} | sediment {r[0]:5.0f} m (Earth {es}) | z {r[1]:+5.0f} m (Earth {ez:+d}) | '
                  f'under the datum {r[2]:.2f} | rock share {r[3]:.2f}')
