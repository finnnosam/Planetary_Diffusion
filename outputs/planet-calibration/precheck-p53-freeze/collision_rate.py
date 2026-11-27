# Pass 54 pre-check: is the thick crust's source (collision) too fast, against Bird's (2003) CCB steps?
import numpy as np, glob, re
from planet_grid.sphere_graph import CubeSphere
from planet_earth.bird import read_dbf
r = read_dbf('data/tectonicplates-master/PB2002_steps.dbf')
ccb = r['STEPCLASS'] == 'CCB'; conv = -r['VELOCITYDI'][ccb]; L = r['STEPLENGTH'][ccb]
def wq(v, w, qs=(25, 50, 75)):
    o = np.argsort(v); c = np.cumsum(w[o]); return [round(float(v[o][np.searchsorted(c, c[-1]*q/100)]), 1) for q in qs]
print('Earth CCB: length %.1fk km, conv p25/50/75 %s mm/yr, length-weighted mean %.1f, shortening area %.3f Mkm2/Myr'
      % (L.sum()/1e3, wq(conv, L), (conv*L).sum()/L.sum(), (np.clip(conv, 0, None)*L).sum()*1e-3/1e6*1e3))
g = CubeSphere(217); h = np.sqrt(g.area.mean())/1e3*0.80     # km of boundary line per node (trench_origin's node factor)
fs = sorted(glob.glob('outputs/planet-plates/freeze-20261031/seed1001/snapshots/t*.npz'), key=lambda f: float(re.search(r't([+-][\d.]+)\.npz', f).group(1)))
rates, lens, area = [], [], []
for f in fs:
    t = float(re.search(r't([+-][\d.]+)\.npz', f).group(1))
    if t < -40: continue
    s = np.load(f); m = s['boundary_type'] == 3
    c = s['convergence_mm_yr'][m].astype(float); rates.append(c); lens.append(m.sum()*h)
    area.append((np.clip(c, 0, None)*h).sum()*1e-3)
c = np.concatenate(rates)
print('planet collision (btype 3), last 40 Myr: length mean %.1fk km, conv p25/50/75 %s mm/yr, mean %.1f, shortening area %.3f Mkm2/Myr'
      % (np.mean(lens)/1e3, wq(c, np.ones_like(c)), c.mean(), np.mean(area)))
print('planet area ledger, collision over 250 Myr: %.3f Mkm2/Myr' % (0.11131947789040604*g.area.sum()/1e12/250))
