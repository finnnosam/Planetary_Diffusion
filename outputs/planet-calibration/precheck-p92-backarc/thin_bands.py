"""Pass 92: the planet's thin crust (total 20-30 km, label-free) by distance from the continent's edge, as
precheck-p91-thin/thin_where.py reads it but without its Earth table (no ETOPO needed: the cloud pre-check). Earth
(CRUST1.0, ice-free; thin_where.txt): 0.81 / 0.68 / 0.37 / 0.17 / 0.098 / 0.058 / 0.009. Args: snapshot tag, seed
folders. Run from the repo root."""
import sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
g = CubeSphere(217); a = g.area
B = [0, 50, 100, 200, 300, 500, 800, 20000]
tag = sys.argv[1]
for sd in sys.argv[2:]:
    q = np.load(f'{sd}/snapshots/{tag}.npz')
    tot = (q['crust_thickness_m']+q['sediment_m'])/1e3; whole = q['continental'].astype(bool) & (tot >= 20)
    d = gl.distance_to(~whole, g.nbr, g.dist)/1e3; thin = whole & (tot < 30)
    z = q['elevation_m']
    inl = whole & (d >= 300); land = inl & (z > 0)
    print(f'{sd} {tag}: thin share {a[thin].sum()/a[whole].sum():.3f}; by band ' +
          ' / '.join(f'{a[thin & (d >= lo) & (d < hi)].sum()/max(a[whole & (d >= lo) & (d < hi)].sum(), 1):.3f}' for lo, hi in zip(B[:-1], B[1:])) +
          f'; mean total by band ' + ' / '.join(f'{np.average(tot[whole & (d >= lo) & (d < hi)], weights=a[whole & (d >= lo) & (d < hi)]):.1f}'
                                               for lo, hi in zip(B[:-1], B[1:]) if (whole & (d >= lo) & (d < hi)).any()) +
          f'; continent area {a[whole].sum()/a.sum():.3f}')
