# Pass 62 diagnostic (no model change): the control's 0 Myr surfaces read with the sea on the datum (Earth's level)
# instead of the planet's low sea: drowned continent, shelf area, land percentiles. Elevations are not re-solved
# (no water reloading); a reading of how much of the shelf / drowned-crust deficit is the sea's level.
import numpy as np, glob
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m
from planet_tectonics.plates import PARAMS
g = CubeSphere(217); a = g.area
def m(T, z):
    cont = T >= 20e3; land = z > 0; sea = ~land
    return [a[cont & sea].sum()/a[cont].sum(), a[sea & (z > -200)].sum()/1e12, a[land].sum()/a.sum(),
            *[np.percentile(z[land], q) for q in (25, 50, 90)]]
print('drowned share | shelf Mkm2 | land frac | land p25 p50 p90')
for name, pat in (('control', 'sweep-20261101-115244/combo-01'), ('combo-02', 'sweep-20261101-191306/combo-02')):
    R0, R1, sl = [], [], []
    for d in sorted(glob.glob(f'/home/claude/pre/sw/{pat}/seed*')):
        zz = np.load(f'{d}/snapshots/t+0000.0.npz'); T = column_thickness_m(zz, PARAMS)
        z = zz['elevation_m'].astype(float); s = float(zz['sea_level_m'])
        # z is relative to the sea; the sea stands s on the datum. On the datum: z + s; sea at the datum: z + s - 0
        R0.append(m(T, z)); R1.append(m(T, z+s)); sl.append(s)
    print(name, 'sea on datum %.0f m' % np.mean(sl))
    print('  as run     ', ' '.join('%8.3f' % x for x in np.mean(R0, 0)))
    print('  sea at 0   ', ' '.join('%8.3f' % x for x in np.mean(R1, 0)))
