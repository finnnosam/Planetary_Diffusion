"""Pass 96 pre-check reader (cloud, seed 1022, planet only): the cover's placement (cover_accommodation.read) at
-250 and 0 Myr, and the drowned band (drowned_band.rows) at 0 Myr, per arm folder."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from cover_accommodation import read, show          # noqa: E402
from drowned_band import rows, BANDS                # noqa: E402
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import sediment_bulk_m
from planet_tectonics.plates import PARAMS

g = CubeSphere(217); a = g.area
for arm in map(Path, sys.argv[1:]):
    sd = arm/'seed1022'
    for t in ('t-0250.0', 't+0000.0'):
        f = np.load(sd/'snapshots'/f'{t}.npz')
        Tc = f['crust_thickness_m'].astype(float); Sb = sediment_bulk_m(f)
        c = (Tc+Sb) >= 20e3
        zd = f['elevation_m'].astype(float)+float(f['sea_level_m'])
        show(f'{arm.name} {t}', read(a[c], Sb[c], Tc[c], zd[c], PARAMS))
        rs = rows(g, Tc, Sb, f['elevation_m'].astype(float))
        print('   continent by elevation about the sea: ' + ' | '.join(
            f"{lo:.0f}..{hi:.0f} {r['share']:.3f} (Tc {r['Tc']:.1f} km, S {r['S']:.0f} m, {r['dist']:.0f} km in)"
            for (lo, hi), r in zip(BANDS, rs)))
