# Pass 60: label-free continental volume (crust + sediment where total >= 20 km), Earth against the planet's start and end
import numpy as np, glob
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, run_params
g = CubeSphere(217); a = g.area
def vol(T): c = T >= 20e3; return (T[c]*a[c]).sum()/1e9/1e6, a[c].sum()/1e12   # Mkm3, Mkm2
e = np.load('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz')
VE, AE = vol(e['T'].astype(float)); print('Earth (CRUST1.0 crsthk, incl. sediment, ice in): %.2f Mkm3 on %.1f Mkm2, mean %.1f km' % (VE, AE, VE/AE))
for d in ['outputs/planet-plates/freeze-20261031/seed1001'] + sorted(glob.glob('outputs/planet-calibration/sweep-20261101-191306/combo-0[24]/seed100[1-6]')):
    p = run_params(d); out = []
    for t in ('-0250.0', '-0040.0', '+0000.0'):
        f = f'{d}/snapshots/t{t}.npz'
        try: V, A = vol(column_thickness_m(np.load(f), p)); out.append('t%s %.2f Mkm3 / %.1f Mkm2 (%.1f km)' % (t, V, A, V/A))
        except FileNotFoundError: pass
    print(d.split('/')[-2] + '/' + d.split('/')[-1] + ': ' + '; '.join(out))
