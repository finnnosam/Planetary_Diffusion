# Pass 63: simultaneously active continental rift length through a run, against Brune et al. (2017, Nat. Geosci. 10,
# 941): > 50,000 km in the Jurassic / Early Cretaceous (Pangea's fragmentation), 20,000 km in the Late Cretaceous, up to
# 30,000 km from the Eocene; East Africa today ~6,000 km of the present network.
# Rift = ridge-boundary nodes (boundary_type 1) on continental crust (the stretched gap crust is continental until
# breakup), read by flag and label-free (crust + sediment >= 20 km). Length = nodes x spacing x f, f from the run's own
# ridge length (stats per step) against its ridge node count, so the double-sided boundary is counted once.
#   python rift_length.py <seed dir> [<seed dir> ...]
import sys, glob, json, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, run_params

g = CubeSphere(217); sp = float(np.sqrt(g.area.mean()))
for s in map(Path, sys.argv[1:]):
    p = run_params(s)
    st = json.load(open(s/'stats.json'))
    ridge_km = {round(float(t)): r for t, r in zip(st.get('t', []), st.get('ridge_km', []))} if 'ridge_km' in st else {}
    print(s, '| t, ridge k km (nodes x sp), continental ridge by flag / label-free (k km; nodes x sp / 2: the boundary is two-sided, checked against the logged ridge length), gap nodes on continental crust (one-sided: nodes x sp)')
    rows = []
    for f in sorted(glob.glob(str(s/'snapshots'/'t*.npz')), key=lambda f: float(Path(f).stem[1:])):
        z = np.load(f); t = float(z['time_myr'])
        if t % 10:
            continue
        rid = z['boundary_type'] == 1
        c = z['continental'].astype(bool)
        lf = column_thickness_m(z, p) >= 20e3
        gap = z['new_crust'].astype(bool)
        r = [t, rid.sum()*sp/2e6, (rid & c).sum()*sp/2e6, (rid & lf).sum()*sp/2e6, (gap & c).sum()*sp/1e6]
        rows.append(r)
        print('%6.0f  %6.1f  %6.1f  %6.1f  %6.1f' % tuple(r))
    R = np.array(rows)
    print('mean over the run: continental ridge by flag %.1f, label-free %.1f k km' % (R[:, 2].mean(), R[:, 3].mean()))
