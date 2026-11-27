# Pass 60: label-free area and thickness bins through a run (every 10 Myr), share of the surface
import sys, glob, re, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, run_params
d = sys.argv[1]
g = CubeSphere(217); A = g.area/g.area.sum()
fs = sorted(glob.glob(f'{d}/snapshots/t*.npz'), key=lambda f: float(re.search(r't([+-][\d.]+)\.npz', f).group(1)))
p = run_params(d)
B = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 99)]
print(' t_myr  label-free  flag   ' + ' '.join('%5d-%-2d' % (lo, min(hi, 99)) for lo, hi in B) + '   flag<20  sea')
for f in fs:
    t = float(re.search(r't([+-][\d.]+)\.npz', f).group(1))
    if t % 10: continue
    s = np.load(f); T = column_thickness_m(s, p); flag = s['continental'].astype(bool)
    print('%6.0f   %.4f   %.4f  ' % (t, A[T >= 20e3].sum(), A[flag].sum()) + ' '.join('%8.4f' % A[(T >= lo*1e3) & (T < hi*1e3)].sum() for lo, hi in B)
          + '   %.4f  %5.0f' % (A[flag & (T < 20e3)].sum(), float(s['sea_level_m'])))
