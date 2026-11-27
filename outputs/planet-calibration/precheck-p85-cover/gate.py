"""Pass 85 (journal 110 §5): the elevation gates pooled over snapshots[-40::5] + the last, as elevation_check reads
them (V2'), and each snapshot's median. Run from the repo root: gate.py <seed folders>"""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import elevation_check as ec
from planet_reviews.crust_review import run_params
g = CubeSphere(217)
for sd in sys.argv[1:]:
    sd = Path(sd); p = run_params(sd)
    snaps = sorted((sd/'snapshots').glob('t*.npz'), key=lambda s: float(s.stem[1:]))
    sel = snaps[-40::5]+[snaps[-1]]
    tests = [ec.snapshot_tests(s, g, p) for s in sel]
    for k in ('young_ocean', 'thin_crust'):
        pooled = ec.weighted_median(np.concatenate([t['_'+k][0] for t in tests]), np.concatenate([t['_'+k][1] for t in tests]))
        print(sd.parent.name, k, 'pooled', round(pooled), 'by snapshot', [round(t[k+'_offset_m']) for t in tests], [s.stem for s in sel][:3], len(sel))
