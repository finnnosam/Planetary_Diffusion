# Pass 62: is the 20-30 km crystalline deficit edge length or width? Label-free continent (column >= 20 km), 46 km.
import numpy as np, glob
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, sediment_bulk_m
from planet_tectonics.plates import PARAMS
g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))
def m(T, S):
    cont = T >= 20e3; X = T-S
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    L = edge.sum()*sp/1e3
    d = gl.distance_to(~cont, nb, g.dist, 1e7)/1e3
    th = cont & (X >= 20e3) & (X < 30e3)
    r = [L/1e3, a[cont].sum()/1e12, a[th].sum()/1e12]
    for lo, hi in ((0, 200), (200, 500), (500, 1e9)):
        r.append(a[th & (d > lo) & (d <= hi)].sum()/1e12)
    r.append(a[th & (d <= 200)].sum()/1e6/L)            # km2 per km of edge within 200 km
    r.append(a[cont & (d <= 200)].sum()/1e12)
    return np.array(r)
hdr = 'edge k km | area | X20-30 | <=200 | 200-500 | >500 km | per km edge (<=200) | cont area <=200 km'
E = np.load('/home/claude/pre/earth_nodes.npz')
print(hdr); print('Earth        ', ' '.join('%7.1f' % x for x in m(E['T'], E['S'])))
for name, pat, t in (('control -250', 'sweep-20261101-115244/combo-01', '-0250.0'), ('control    0', 'sweep-20261101-115244/combo-01', '+0000.0'), ('combo-02   0', 'sweep-20261101-191306/combo-02', '+0000.0')):
    R = []
    for d in sorted(glob.glob(f'/home/claude/pre/sw/{pat}/seed*')):
        z = np.load(f'{d}/snapshots/t{t}.npz')
        R.append(m(column_thickness_m(z, PARAMS), sediment_bulk_m(z)))
    print(name, ' '.join('%7.1f' % x for x in np.mean(R, 0)))
