"""Pass 94 pre-check reader (seed 1007, control vs 'column', 0 Myr and the -40...0 window): land fraction, continents
(crust + bulk sediment >= 30 km, continental flag) on the datum p25 / p50 / p75 and share under it, drowned flagged
continent, the -1,000 ... -200 m band about the sea, land p50, the sea on the datum."""
import glob, sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.isostasy import sediment_column
from planet_reviews.crust_review import run_params
g = CubeSphere(217); A = g.area/g.area.sum()
def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); c = (c-.5*w[o])/c[-1]; return np.interp(q, c, x[o])
for arm in sys.argv[1:]:
    run = f'{arm}/seed1007'; p = run_params(run); rows = []
    fs = sorted(glob.glob(f'{run}/snapshots/t*.npz'), key=lambda s: float(s.split('/')[-1][1:-4]))
    for f in [x for x in fs if float(x.split('/')[-1][1:-4]) >= -40]:
        d = np.load(f); sl = float(d['sea_level_m']); z = d['elevation_m'].astype(float); zd = z+sl
        cont = d['continental'].astype(bool); sh = sediment_column(d['sediment_m'].astype(float), cont, p)[0]
        tot = (d['crust_thickness_m']+sh)/1e3; c30 = cont & (tot >= 30)
        land = z > 0
        rows.append([A[land].sum(), *wq(zd[c30], A[c30], [.25, .5, .75]), A[c30 & (zd < 0)].sum()/A[c30].sum(),
                     A[cont & ~land].sum()/A[cont].sum(), A[(z >= -1000) & (z < -200)].sum(), wq(z[land], A[land], .5), sl])
    r = np.array(rows)
    for nm, v in (('0 Myr', r[-1]), ('window mean', r.mean(0))):
        print(f'{arm.split("/")[-1]:8s} {nm:11s}: land {v[0]:.3f} | continents on the datum p25/50/75 {v[1]:.0f} / {v[2]:.0f} / {v[3]:.0f}, under {v[4]:.3f} | drowned flagged {v[5]:.2f} | -1000..-200 {v[6]:.3f} | land p50 {v[7]:.0f} | sea {v[8]:+.0f}')
