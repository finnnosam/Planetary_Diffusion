"""Pass 87 (as pass 85's gate_cause.py, journal 110 §5): a snapshot's young-ocean gate offset
(elevation_check.snapshot_tests) by quartile of the sediment within ~140 km (three rings of neighbours), and by
quartile of the crust thickness within the same rings (docked arcs thicken margins). Run from the repo root:
gate_cause.py <seed folder> <snapshot stem ...>"""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import elevation_check as ec
from planet_reviews.crust_review import run_params
from planet_landscape import graph_landscape as gl
import planet_tectonics.plates as P
g = CubeSphere(217); a = g.area
sd = Path(sys.argv[1]); p = run_params(sd)
for stem in sys.argv[2:]:
    f = sd/'snapshots'/f'{stem}.npz'; z = np.load(f)
    t = ec.snapshot_tests(f, g, p)
    off, w = t['_young_ocean']
    cont = z['continental'].astype(bool); S = z['sediment_bulk_m'].astype(float)
    age = z['ocean_age_myr']; T = z['crust_thickness_m'].astype(float)
    d = gl.distance_to(cont, g.nbr, g.dist, 150e3)
    sedh = P.sediment_column(z['sediment_m'].astype(float), cont, p)[0]
    yo = (~cont) & (age < 20) & (T-p['ocean_thick_m'] < 1000) & (sedh < 500) & (d <= 100e3)
    assert yo.sum() == len(off), (yo.sum(), len(off))
    nb = g.nbr; idx = np.where(yo)[0]
    near, thick = [], []
    for i in idx:
        ring = set([i]); fr = {i}
        for _ in range(3):
            fr = {j for k in fr for j in nb[k] if j >= 0} - ring; ring |= fr
        ring = np.array(sorted(ring)); near.append((S[ring]*a[ring]).sum()/a[ring].sum()); thick.append(((T+S)[ring]*a[ring]).sum()/a[ring].sum())
    near, thick = np.array(near), np.array(thick)
    q = np.array_split(np.argsort(near), 4); qt = np.array_split(np.argsort(thick), 4)
    print(sd.parent.name, stem, 'young nodes', yo.sum(), 'offset median', round(float(np.median(off))),
          '| by sediment within ~140 km: ' + '; '.join(f'{near[k].mean():.0f} m -> {np.median(off[k]):.0f}' for k in q)
          + ' | by column within ~140 km: ' + '; '.join(f'{thick[k].mean()/1e3:.1f} km -> {np.median(off[k]):.0f}' for k in qt))
