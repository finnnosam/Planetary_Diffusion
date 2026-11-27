# Pass 64: how wide is the thickened crust beside Earth's collisions? Bird (2003) CCB steps (continental convergent
# boundaries) and orogen polygons; CRUST1.0 total thickness on n = 217. Thick = >= 45 km inside Bird's orogens and
# within 2,000 km of a CCB step; width = its area per km of CCB, and the distance of its nodes from the CCB.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_grid.paths import DATA
from planet_landscape import graph_landscape as gl
from planet_earth.bird import read_dbf, xyz_from_lonlat, orogen_mask
from planet_reviews.margin_width import earth_nodes
g = CubeSphere(217); a = g.area
s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
ccb = s['STEPCLASS'] == 'CCB'
L = s['STEPLENGTH'][ccb].sum()
pts = []
for lo0, la0, lo1, la1 in zip(s['STARTLONG'][ccb], s['STARTLAT'][ccb], s['FINALLONG'][ccb], s['FINALLAT'][ccb]):
    A, B = xyz_from_lonlat(lo0, la0), xyz_from_lonlat(lo1, la1)
    for t in np.linspace(0, 1, 6):
        p = A*(1-t)+B*t; pts.append(p/np.linalg.norm(p))
front = np.zeros(g.size, bool); front[fields.nearest_nodes(g, np.array(pts))] = True
T, _ = earth_nodes(g)
oro = orogen_mask(g, DATA/'tectonicplates-master'/'GeoJSON'/'PB2002_orogens.json')
d = gl.distance_to(front, g.nbr, g.dist, 3000e3)/1e3
print(f'CCB length {L:,.0f} km (Bird steps), front nodes {front.sum()}')
for lim in (45e3, 50e3, 55e3):
    for name, sel in (('in orogens', oro), ('anywhere', np.ones(g.size, bool))):
        th = (T >= lim) & sel & (d <= 2000)
        q = np.percentile(d[th], [25, 50, 75, 90])
        print(f'crust >= {lim/1e3:.0f} km, {name}, <= 2,000 km of a CCB: {a[th].sum()/1e12:5.2f} Mkm2 = {a[th].sum()/1e6/L:5.0f} km^2 per km of CCB; '
              f'distance from the CCB p25/50/75/90 {q[0]:.0f} / {q[1]:.0f} / {q[2]:.0f} / {q[3]:.0f} km')
# the kernel: share of the continental crust (total >= 20 km) inside Bird's orogens that is >= 45 km thick, by distance
# from the CCB (100 km bins to 2,000 km); the collision spread's weight per unit area (journal 86)
cont = T >= 20e3
edges = np.arange(0, 2001, 100)
sh = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = cont & oro & (d >= lo) & (d < hi)
    sh.append(a[m & (T >= 45e3)].sum()/max(a[m].sum(), 1))
print('kernel (share >= 45 km by 100 km bin from the CCB, 0-2,000 km):', ', '.join('%.3f' % x for x in sh))
