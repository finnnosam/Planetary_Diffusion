# Pass 64: where Earth's crust thickens beside its collisions today. GSRM v2.1 (Kreemer et al. 2014) areal dilatation
# exx + eyy < 0 on continental crust (CRUST1.0 total >= 20 km) is the thickening rate (incompressible crust:
# dT/dt = -T (exx + eyy)); its volume rate by distance from Bird's (2003) CCB steps, leaving out crust within 300 km
# of a SUB step (locked megathrusts' interseismic elastic shortening, journal 82). Gives the collision spread's kernel.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_grid.paths import DATA
from planet_landscape import graph_landscape as gl
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_earth.deformation import gsrm_nodes
from planet_reviews.margin_width import earth_nodes
g = CubeSphere(217); a = g.area
s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
def front(cls):
    k = s['STEPCLASS'] == cls; pts = []
    for lo0, la0, lo1, la1 in zip(s['STARTLONG'][k], s['STARTLAT'][k], s['FINALLONG'][k], s['FINALLAT'][k]):
        A, B = xyz_from_lonlat(lo0, la0), xyz_from_lonlat(lo1, la1)
        for t in np.linspace(0, 1, 6):
            p = A*(1-t)+B*t; pts.append(p/np.linalg.norm(p))
    m = np.zeros(g.size, bool); m[fields.nearest_nodes(g, np.array(pts))] = True
    return m, s['STEPLENGTH'][k].sum(), (s['STEPLENGTH'][k]*np.abs(s['VELOCITYDI'][k] if 'VELOCITYDI' in s else 0)).sum()
ccb, L, _ = front('CCB'); sub, _, _ = front('SUB')
T, _ = earth_nodes(g); cont = T >= 20e3
exx, eyy, exy, cover = gsrm_nodes(g)
rate = np.clip(-(exx+eyy), 0, None)*T*a           # m^3/yr of thickening at each node
d = gl.distance_to(ccb, g.nbr, g.dist, 4000e3)/1e3
dsub = gl.distance_to(sub, g.nbr, g.dist, 4000e3)/1e3
sel = cont & (dsub > 300) & (d <= 3000)
print(f'CCB {L:,.0f} km; thickening on continental crust > 300 km from SUB and <= 3,000 km of a CCB: {rate[sel].sum()/1e9:.3f} km^3/yr '
      f'({rate[sel].sum()/1e9/L*1e3:.1f} km^3/Myr per km of CCB)')
edges = np.array([0, 100, 200, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000, 3000])
tot = rate[sel].sum(); cum = 0.
print('distance from the CCB (km) | share of the thickening | cumulative | per unit continental area (relative to 0-100 km)')
dens0 = None
for lo, hi in zip(edges[:-1], edges[1:]):
    m = sel & (d >= lo) & (d < hi)
    sh = rate[m].sum()/tot; cum += sh
    dens = rate[m].sum()/max(a[m].sum(), 1); dens0 = dens0 or dens
    print(f'{lo:5.0f}-{hi:<5.0f} {sh:6.3f} {cum:6.3f} {dens/dens0:6.3f}')
w = rate[sel]; dd = d[sel]; o = np.argsort(dd); c = np.cumsum(w[o])/w.sum()
print('distance quantiles of the thickening p25/50/75/90: ' + ' / '.join('%.0f' % dd[o][np.searchsorted(c, q)] for q in (.25, .5, .75, .9)) + ' km')
