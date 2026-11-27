import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.elevation_check import snapshot_tests, weighted_median
from planet_reviews.crust_review import run_params
g = CubeSphere(217)
for arm in sys.argv[1:]:
    sd = Path(arm)/'seed1020'; p = run_params(sd)
    snaps = sorted((sd/'snapshots').glob('t*.npz'), key=lambda s: float(s.stem[1:]))
    tests = [snapshot_tests(s, g, p) for s in snaps[-40::5]+[snaps[-1]]]
    o = {k: weighted_median(np.concatenate([t['_'+k][0] for t in tests]), np.concatenate([t['_'+k][1] for t in tests])) for k in ('young_ocean', 'thin_crust')}
    print(f"{sd}: gates pooled young-ocean {o['young_ocean']:+.0f} m, thin-crust {o['thin_crust']:+.0f} m; young ocean above sea max {max(t['young_ocean_above_sea'] for t in tests):.3f}")
