# Pass 62: edge history on the freeze run (seed 1001): label-free edge, passive edge (> 300 km from any boundary),
# crystalline 20-30 km crust near passive edges, through the run.
import numpy as np, glob, json
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, sediment_bulk_m
from planet_tectonics.plates import PARAMS
g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))
def m(T, S, bnd):
    cont = T >= 20e3; X = T-S
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    far = gl.distance_to(bnd, nb, g.dist, 5000e3) > 300e3
    pas = edge & far
    dpe = gl.distance_to(pas, nb, g.dist, 600e3)
    thin = cont & (X < 30e3)
    return [edge.sum()*sp/1e6, pas.sum()*sp/1e6, (edge & ~far).sum()*sp/1e6, a[cont].sum()/1e12, a[thin].sum()/1e12,
            a[thin & (dpe <= 500e3)].sum()/1e12, a[thin & (dpe <= 500e3)].sum()/max(pas.sum()*sp, 1)/1e3]
print('t     edge  passive  active(k km)  cont   X<30   X<30 near passive  per km passive')
for f in sorted(glob.glob('/home/claude/pre/fz/snapshots/*.npz'), key=lambda f: float(f.split('/t')[-1][:-4])):
    z = np.load(f)
    r = m(column_thickness_m(z, PARAMS), sediment_bulk_m(z), z['boundary_type'] != 0)
    print('%6.0f ' % z['time_myr'] + ' '.join('%7.1f' % x for x in r))
E = np.load('/home/claude/pre/earth_nodes.npz')
from planet_earth.bird import earth_boundaries
print('Earth  ' + ' '.join('%7.1f' % x for x in m(E['T'], E['S'], earth_boundaries(g))))
ev = json.load(open('/home/claude/pre/fz/stats.json'))['event_list']
for e in ev: print(e)
