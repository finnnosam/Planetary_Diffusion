# Pass 62: label-free continental pieces by size: count, area, edge, passive edge (> 300 km from boundaries); Earth vs
# the control at 0 Myr (combo-01 of sweep-20261101-115244, six seeds) and combo-02 of sweep-20261101-191306.
import numpy as np, glob
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics.graphops import components
from planet_reviews.crust_review import column_thickness_m
from planet_tectonics.plates import PARAMS
from planet_earth.bird import earth_boundaries
g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))
CL = [(0, .1), (.1, 1), (1, 10), (10, 1e9)]
def m(T, bnd):
    cont = T >= 20e3
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    far = gl.distance_to(bnd, nb, g.dist, 5000e3) > 300e3
    lab = components(cont, nb)
    A = np.bincount(lab[cont], weights=a[cont])/1e12
    E = np.bincount(lab[edge], minlength=len(A))*sp/1e6
    P = np.bincount(lab[edge & far], minlength=len(A))*sp/1e6
    out = []
    for lo, hi in CL:
        s = (A >= lo) & (A < hi)
        out += [s.sum(), A[s].sum(), E[s].sum(), P[s].sum()]
    return np.array(out)
print('class Mkm2: ' + ' | '.join(f'{lo}-{hi}: n area edge passive(k km)' for lo, hi in CL))
E = np.load('/home/claude/pre/earth_nodes.npz')
print('Earth   ', ' '.join('%6.1f' % x for x in m(E['T'], earth_boundaries(g))))
for name, pat in (('control', 'sweep-20261101-115244/combo-01'), ('combo-02', 'sweep-20261101-191306/combo-02')):
    R = [m(column_thickness_m(z := np.load(f'{d}/snapshots/t+0000.0.npz'), PARAMS), z['boundary_type'] != 0) for d in sorted(glob.glob(f'/home/claude/pre/sw/{pat}/seed*'))]
    print('%-8s' % name, ' '.join('%6.1f' % x for x in np.mean(R, 0)))
