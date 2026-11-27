"""Pass 97 probe: the setting of the drowned continent (label-free continent, crust + bulk sediment >= 20 km, at
-1000 ... -200 m about its sea) at 0 Myr against the rest of the continent: share within 600 / 1,000 km of a subduction
boundary node, lithosphere age, crystalline crust, cumulative tectonic thickening (snapshot field), and the share whose
crust is thinner than 32 km. Six-seed area means. Read, not a rule. Usage (repo root): python3 drowned_setting.py <arm> ..."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import sediment_bulk_m
g = CubeSphere(217); a = g.area; SUB = 2
for arm in sys.argv[1:]:
    R = []
    for sd in sorted(Path(arm).glob('seed*')):
        q = np.load(sd/'snapshots'/'t+0000.0.npz')
        Tc = q['crust_thickness_m']; S = q['sediment_bulk_m']; z = q['elevation_m']
        c = (Tc+S) >= 20e3
        ds = gl.distance_to(q['boundary_type'] == SUB, g.nbr, g.dist)/1e3
        L = q['lithosphere_age_myr']; th = q['tectonic_thickening_m']
        row = []
        for m in (c & (z >= -1000) & (z < -200), c & (z >= 0)):
            w = a[m]; W = w.sum()
            row += [W/a.sum(), (w*(ds[m] < 600)).sum()/W, (w*(ds[m] < 1000)).sum()/W, np.nansum(w*L[m])/W,
                    (w*Tc[m]).sum()/W/1e3, (w*th[m]).sum()/W/1e3, (w*(Tc[m] < 32e3)).sum()/W]
        R.append(row)
    m = np.mean(R, 0)
    for k, lab in ((0, 'drowned -1000..-200'), (7, 'land')):
        print(f'{Path(arm).name} {lab}: share {m[k]:.3f}; within 600 / 1000 km of subduction {m[k+1]:.2f} / {m[k+2]:.2f}; '
              f'L age {m[k+3]:.0f} Myr; Tc {m[k+4]:.1f} km; tectonic thickening {m[k+5]:+.2f} km; Tc < 32 km {m[k+6]:.2f}')
