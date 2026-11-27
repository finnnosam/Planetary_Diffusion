"""Pass 97 probe: the submerged continent's depth by distance from land. Label-free continent (crust + bulk sediment
>= 20 km) under the sea; distance (km, graph) from the nearest land node (elevation >= 0 about the sea); area-weighted
p25 / p50 / p75 depth and the share of the submerged continent per distance bin. Earth: ETOPO on the grid, CRUST1.0,
ice-sheet nodes out of the sample (land under ice still counts as land). Planet: 0 Myr snapshots, six-seed means.
Also the same for the label-free continent beyond 300 km from any plate boundary (Earth: Bird PB2002 nodes from
earth-plates-n217.npz if present). Read, not a rule. Usage (repo root): python3 shelf_profile.py <arm> ..."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.coast_budget import wmed
from planet_reviews.crust_review import sediment_bulk_m
BINS = (0, 50, 100, 150, 200, 300, 400, 600, 1000, 3000)
g = CubeSphere(217); a = g.area


def wq(v, w, q):
    o = np.argsort(v); c = np.cumsum(w[o])
    return v[o][np.searchsorted(c, q*c[-1])] if c.size and c[-1] > 0 else np.nan


def prof(cont, z, keep):
    sub = cont & (z < 0) & keep
    d = gl.distance_to(z >= 0, g.nbr, g.dist)/1e3
    out = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        m = sub & (d >= lo) & (d < hi)
        out.append([a[m].sum()/a[sub].sum(), wq(z[m], a[m], .25), wq(z[m], a[m], .5), wq(z[m], a[m], .75)])
    return np.array(out), a[sub].sum()/a.sum()


def show(name, P, share):
    print(f'{name}: submerged continent {share:.3f} of the surface')
    print('   ' + ' | '.join(f'{lo}-{hi} km: {r[0]:.2f} {r[1]:.0f}/{r[2]:.0f}/{r[3]:.0f}' for (lo, hi), r in
                          zip(zip(BINS[:-1], BINS[1:]), P)))


lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
from planet_earth.reference import ice_mask
z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
ice = ice_mask(g)
show('Earth (share per bin, depth p25/p50/p75 m)', *prof(T >= 20e3, z, ~ice))
for arm in map(Path, sys.argv[1:]):
    acc, sh = [], []
    for s in sorted(arm.glob('seed*')):
        f = np.load(s/'snapshots'/'t+0000.0.npz')
        cont = (f['crust_thickness_m']+sediment_bulk_m(f)) >= 20e3
        P, x = prof(cont, f['elevation_m'].astype(float), np.ones(g.size, bool))
        acc.append(P); sh.append(x)
    show(f'{arm.name} (six-seed means)', np.nanmean(acc, 0), float(np.mean(sh)))
print('planet graded shelf ceiling today: -(10 + 1.2 m/km d) to 150 km from land: 50 km -70, 100 km -130, 150 km -190; none beyond')
