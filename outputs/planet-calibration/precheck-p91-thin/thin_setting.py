"""Pass 91: Earth's thin continental crust (CRUST1.0 total 20-30 km, non-oceanic, ice-free) by distance from the
continent's edge and by distance to Bird 2003's nearest subduction step (overriding or not) and to any boundary.
Run from the repo root."""
import numpy as np
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import lithosphere as lt
from planet_earth.reference import ice_mask
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_grid import fields
g = CubeSphere(217)
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]
tn = np.array([l.split() for l in (DATA/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
whole = ~np.isin(tn, lt.OCEANIC_TYPES) & (T >= 20)
keep = whole & ~ice_mask(g)
d = gl.distance_to(~whole, g.nbr, g.dist)/1e3
s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
print('step classes', sorted(set(s['STEPCLASS'])))
cls = s['STEPCLASS']
mid = xyz_from_lonlat(.5*(s['STARTLONG']+s['FINALLONG']), .5*(s['STARTLAT']+s['FINALLAT']))
def dist_to(mask_steps):
    m = np.zeros(g.size, bool); m[fields.nearest_nodes(g, mid[mask_steps])] = True
    return gl.distance_to(m, g.nbr, g.dist)/1e3
dsub = dist_to(cls == 'SUB'); dany = dist_to(np.ones(len(cls), bool))
thin = keep & (T < 30); a = g.area
for lo, hi in ((0, 200), (200, 20000)):
    k = thin & (d >= lo) & (d < hi); A = a[k].sum()
    print(f'\nthin crust {lo}-{hi} km from the edge: {A/a[thin].sum():.2f} of all thin')
    for name, dd in (('subduction step', dsub), ('any boundary', dany)):
        print('  '+name+': ' + ', '.join(f'<{x} km {a[k & (dd < x)].sum()/A:.2f}' for x in (300, 600, 1000, 1500)))
    nk = keep & (d >= lo) & (d < hi); An = a[nk].sum()
    print('  (all continent in this band: subduction <600 km %.2f, <1000 km %.2f)' % (a[nk & (dsub < 600)].sum()/An, a[nk & (dsub < 1000)].sum()/An))
