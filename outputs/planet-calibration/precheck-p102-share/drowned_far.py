"""Pass 102 probe: the drowned normal crust (label-free continent, crystalline >= 30 km, -1000 ... -200 m about its
sea) split by distance to subduction (< 1000 / >= 1000 km), at 0 Myr: share of the surface, area-mean dyn (air-loaded
term as continent.py), zs and zs - dyn, and the share of it that would stand above -200 m without the field. Beside:
crystalline >= 30 km continent under its sea with and without the field. Read, not a rule.
Usage (repo root): python3 drowned_far.py <arm> ..."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'precheck-p82-stepback'))
import common as C
from planet_landscape import graph_landscape as gl
for arm in sys.argv[1:]:
    R = []
    for s in C.SEEDS:
        P = C.planet(s, 't+0000.0', arm); g = P['g']; a = P['a']
        q = np.load(Path(arm)/f'seed{s}'/'snapshots'/'t+0000.0.npz')
        Tc = q['crust_thickness_m'].astype(float); Sb = q['sediment_bulk_m'].astype(float)
        c = ((Tc+Sb) >= 20e3) & (Tc >= 30e3); zs = P['zd']-P['sea']; dyn = P['dyn']
        ds = gl.distance_to(q['boundary_type'] == 2, g.nbr, g.dist)/1e3
        dr = c & (zs >= -1000) & (zs < -200)
        row = []
        for m in (dr & (ds < 1000), dr & (ds >= 1000)):
            w = a[m]; W = max(w.sum(), 1e-30)
            row += [W/a.sum(), (w*dyn[m]).sum()/W, (w*zs[m]).sum()/W, (w*(zs-dyn)[m]).sum()/W,
                    (w*((zs-dyn)[m] >= -200)).sum()/W]
        wc = a[c]; row += [(wc*(zs[c] < 0)).sum()/wc.sum(), (wc*((zs-dyn)[c] < 0)).sum()/wc.sum()]
        R.append(row)
    M = np.mean(R, 0)
    print(Path(arm).name, '(6 seeds, 0 Myr, area means)')
    for i, lab in enumerate(('< 1000 km', '>= 1000 km')):
        sh, dy, z, zd, up = M[5*i:5*i+5]
        print(f'  drowned Tc>=30, {lab:10} from subduction: share {sh:.4f} | dyn {dy:+5.0f} | zs {z:+5.0f} | zs-dyn {zd:+5.0f} m'
              f' | above -200 m without the field {up:.2f}')
    print(f'  crystalline >= 30 km continent under its sea: with the field {M[10]:.3f}, without {M[11]:.3f}')
