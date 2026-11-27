"""Pass 98 probe: the drowned continent (label-free, crust + bulk sediment >= 20 km, -1000 ... -200 m about its sea)
split by crystalline crust < 30 / >= 30 km, at 0 Myr. Per part: share of the surface, area means of the elevation
law's terms (continent.py's: airy, dens, sedc, C, dyn, rest) and of the crust, sediment, lithosphere thickness, the
median distance from the ocean and the share within 600 km of subduction. Land (>= 0 m) beside it. Read, not a rule.
Usage (repo root): python3 drowned_split.py <arm> ..."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'precheck-p82-stepback'))
import common as C
from planet_reviews.coast_profile import distance_to_ocean_km
from planet_landscape import graph_landscape as gl
K = ['share', 'Tc', 'S', 'L', 'dist', 'sub600', 'zs', 'airy', 'dens', 'sedc', 'C', 'dyn', 'rest']
for arm in sys.argv[1:]:
    R = []
    for s in C.SEEDS:
        P = C.planet(s, 't+0000.0', arm); g = P['g']; a = P['a']
        q = np.load(Path(arm)/f'seed{s}'/'snapshots'/'t+0000.0.npz')
        Tc = q['crust_thickness_m'].astype(float); Sb = q['sediment_bulk_m'].astype(float)
        c = (Tc+Sb) >= 20e3; zs = P['zd']-P['sea']
        d = distance_to_ocean_km(g, c); ds = gl.distance_to(q['boundary_type'] == 2, g.nbr, g.dist)/1e3
        dr = c & (zs >= -1000) & (zs < -200)
        row = []
        for m in (dr & (Tc < 30e3), dr & (Tc >= 30e3), c & (zs >= 0)):
            w = a[m]; W = w.sum(); mean = lambda x: float((w*x[m]).sum()/W)
            o = np.argsort(d[m]); cw = np.cumsum(w[o])
            row.append([W/a.sum(), mean(Tc)/1e3, mean(Sb), mean(P['L']), d[m][o][np.searchsorted(cw, .5*cw[-1])],
                        float((w*(ds[m] < 600)).sum()/W), mean(zs), mean(P['airy']), mean(P['dens']), mean(P['sedc']),
                        mean(P['cterm']), mean(P['dyn']), mean(P['rest'])])
        R.append(row)
    M = np.mean(R, 0)
    print(Path(arm).name, '(6 seeds, 0 Myr, area means; terms on the datum as continent.py; zs about the sea)')
    print('  part              | ' + ' | '.join(f'{k:>6}' for k in K))
    for lab, r in zip(('drowned Tc < 30', 'drowned Tc >= 30', 'land'), M):
        print(f'  {lab:17} | ' + ' | '.join(f'{v:6.3f}' if k in ('share', 'sub600') else f'{v:6.1f}' if k == 'Tc' else f'{v:6.0f}' for k, v in zip(K, r)))
