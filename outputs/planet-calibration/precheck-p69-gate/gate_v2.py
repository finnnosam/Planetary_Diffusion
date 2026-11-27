"""Pass 69: the elevation gates (V2') of a sweep with the corrected checker (C's term in the column's own level)."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import elevation_check as ec
from planet_reviews.crust_review import run_params

g = CubeSphere(217)
for arm in sorted(Path(sys.argv[1]).glob('*/')):
    for sd in sorted(arm.glob('seed*')):
        p = run_params(sd)
        snaps = sorted((sd/'snapshots').glob('t*.npz'), key=lambda s: float(s.stem[1:]))
        t = [ec.snapshot_tests(s, g, p) for s in snaps[-40::5]+[snaps[-1]]]
        r = {k: ec.weighted_median(np.concatenate([x['_'+k][0] for x in t]), np.concatenate([x['_'+k][1] for x in t]))
             for k in ('young_ocean', 'thin_crust')}
        above = max(x['young_ocean_above_sea'] for x in t)
        ok = all(abs(v) <= 300 for v in r.values()) and np.mean([x['young_ocean_above_sea'] for x in t]) <= .01
        print(f"{arm.name}/{sd.name}: young-ocean {r['young_ocean']:+.0f} m, thin-crust {r['thin_crust']:+.0f} m, "
              f"young ocean above sea max {above:.3f} -> {'pass' if ok else 'FAIL'}", flush=True)
