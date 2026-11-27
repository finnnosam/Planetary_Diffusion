"""Pass 89 (journal 115): T3's continental-side trench split by the origin of its continental side. Earth's T3 (Bird's
SUB steps, CRUST1.0 types 100 and 200 km to both sides; planet_reviews.trench_origin) counts island arc (J-) and
forearc (K-) crust as oceanic, so arcs are intra-oceanic there; the planet's T3 counts any continental-flag node, so arc
crust turned continental (>= 20 km) types its trench as continental-side (journal 114 sec. 4: likely T3's rise since D).
For each subduction contour segment of T3 (planet_reviews.trench_origin.snapshot_metrics' cs), the continental side's
crust: made during the run (crust age < time since the start) or the start's, and the continental piece it lies on
(connected continental nodes): < 1 Mkm2, 1-10, >= 10. Window means over -100 ... 0 Myr every 10 Myr (T3's window),
per arm (six seeds pooled as seed means). Read only.
Run from the repo root: t3_origin.py <arm folder> [...]"""
import sys, json
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.graphops import components
from planet_reviews import boundary_lengths as bl
from planet_reviews import trench_origin as to
g = CubeSphere(217)
WINDOW = tuple(float(t) for t in range(-100, 1, 10))
def seeds_of(d):
    d = Path(d); return [d] if (d/'snapshots').exists() else sorted(d.glob('seed*'))
def snap(sd, t):
    f = to.snap_path(sd, t)
    if not Path(str(f)+'.npz').exists():
        return None
    z = np.load(str(f)+'.npz'); plates = json.loads(Path(str(f)+'-plates.json').read_text())
    pid = z['plate_id']; cont = z['continental'].astype(bool)
    s = bl.typed_segments(g, pid, plates, cont)
    L, typ, na, nb = s['L'], s['typ'], s['na'], s['nb']
    bts = z['boundary_type']
    side_sub = (bts[na] == P.SUBDUCTION) | (bts[nb] == P.SUBDUCTION)
    noslab = (typ == 2) & ((bts[na] == P.CONVERGENT) | (bts[nb] == P.CONVERGENT)) & ~side_sub
    oo = ~(cont[na] | cont[nb])
    cs = (typ == 2) & ~oo & ~noslab
    side = np.where(cont[na], na, nb)                    # the continental side (either if both)
    age = z['crust_age_myr'].astype(float)
    since = float(z['time_myr'])+250.
    made = age[side] < since
    lab = components(cont, g.nbr)
    sz = np.bincount(lab[lab >= 0], weights=g.area[lab >= 0])/1e12
    psz = sz[np.clip(lab[side], 0, None)]
    th = z['crust_thickness_m'].astype(float)[side]
    r = {'T3': L[cs].sum(), 'made': L[cs & made].sum(), 'start': L[cs & ~made].sum(),
         'piece<1': L[cs & (psz < 1)].sum(), 'piece1-10': L[cs & (psz >= 1) & (psz < 10)].sum(),
         'piece>=10': L[cs & (psz >= 10)].sum(), 'made_piece<10': L[cs & made & (psz < 10)].sum(),
         'made_thin<30km': L[cs & made & (th < 30e3)].sum()}
    return r
for arm in sys.argv[1:]:
    rows = []
    for sd in seeds_of(arm):
        w = [r for r in (snap(sd, t) for t in WINDOW) if r]
        if w:
            rows.append({k: np.mean([r[k] for r in w])/1e3 for k in w[0]})
    if not rows:
        print(f'{arm}: no snapshots'); continue
    m = {k: np.mean([r[k] for r in rows]) for k in rows[0]}
    print(f"{Path(arm).name} ({len(rows)} seeds; window -100..0, k km): T3 {m['T3']:.1f} = continental side made in the run "
          f"{m['made']:.1f} + the start's {m['start']:.1f}; by piece < 1 / 1-10 / >= 10 Mkm2 {m['piece<1']:.1f} / "
          f"{m['piece1-10']:.1f} / {m['piece>=10']:.1f}; made in the run on pieces < 10 Mkm2 {m['made_piece<10']:.1f}, "
          f"made with crust < 30 km {m['made_thin<30km']:.1f} (Earth 30.7k, island arcs counted intra-oceanic)")
