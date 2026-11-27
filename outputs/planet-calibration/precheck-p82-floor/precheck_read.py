"""Pass 82 (journal 107): read the cloud pre-check (seed 1013; today's PARAMS and ocean_depth_rule 'earth'). The runs
were stopped between -10 and -5 Myr of 250 (the workspace paused), so the read is at -10 Myr and no stats.json exists: the
parameters are PARAMS + the arm's rule. Land share and p25 / p50 / p90 (area-weighted), sea on the datum, room below
the datum, the highest node over the saves, the elevation gates on the last snapshot. Read only, not pooled.
Run from the repo root: precheck_read.py <runs folder>"""
import sys, glob, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_reviews import elevation_check as ec
g = CubeSphere(217); a = g.area/g.area.sum()
def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); c = (c-.5*w[o])/c[-1]; return np.interp(q, c, x[o])
for name, rule in (('control', 'gdh1'), ('earth', 'earth')):
    sd = f'{sys.argv[1]}/{name}/seed1013/snapshots'; p = {**P.PARAMS, 'ocean_depth_rule': rule}
    zmax = max(float(np.load(f)['elevation_m'].max()) for f in sorted(glob.glob(sd+'/t*.npz')))
    for tag in ('t-0250.0', 't-0100.0', 't-0010.0'):
        q = np.load(f'{sd}/{tag}.npz'); z = q['elevation_m'].astype(float); sl = float(q['sea_level_m']); zd = z+sl; land = z > 0
        T = (q['crust_thickness_m']+q['sediment_m'])/1e3; c20 = T >= 20; fl = (T < 10) & (zd < 0)
        d = (T >= 30) & q['continental'].astype(bool)
        print(f"{name} {tag}: sea {sl:.0f} m; land {a[land].sum():.3f}, p25/p50/p90 " + '/'.join(f'{v:.0f}' for v in wq(z[land], a[land], [.25, .5, .9]))
              + f"; room below the datum {np.sum(a*np.maximum(-zd, 0)):.0f} m; floor (crust < 10 km) share {a[fl].sum():.3f} mean depth {np.sum(a[fl]*-zd[fl])/a[fl].sum():.0f} m; "
              f"continent >= 30 km on the datum p25/p50/p75 " + '/'.join(f'{v:.0f}' for v in wq(zd[d], a[d], [.25, .5, .75])) + f"; crust > 45 km {a[T > 45].sum()/a[c20].sum():.3f}; label-free area {a[c20].sum():.3f}")
    t = ec.snapshot_tests(f'{sd}/t-0010.0.npz', g, p)
    print(f"{name}: highest node over the saves {zmax:.0f} m; gates at -10 Myr: " + ', '.join(f'{k} {v:.3f}' if abs(v) < 5 else f'{k} {v:.0f}' for k, v in t.items() if isinstance(v, float)))
