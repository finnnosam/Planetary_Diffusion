"""Pass 92 (c): Earth's thin continental crust (CRUST1.0 total 20-30 km, non-oceanic, ice-free) by distance to Bird 2003's
nearest subduction step: thin share per 100 km bin against the continent's share far from subduction, and where the
excess thin area lies (the back-arc band's distance and width). Run from the repo root."""
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
d_edge = gl.distance_to(~whole, g.nbr, g.dist)/1e3
s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
cls = s['STEPCLASS']
mid = xyz_from_lonlat(.5*(s['STARTLONG']+s['FINALLONG']), .5*(s['STARTLAT']+s['FINALLAT']))
m = np.zeros(g.size, bool); m[fields.nearest_nodes(g, mid[cls == 'SUB'])] = True
dsub = gl.distance_to(m, g.nbr, g.dist)/1e3
a = g.area; thin = keep & (T < 30)
for lab, sel in (('all continent', keep), ('beyond 200 km of the edge', keep & (d_edge >= 200)),
                 ('within 200 km of the edge', keep & (d_edge < 200))):
    far = sel & (dsub >= 1500); base = a[far & thin].sum()/max(a[far].sum(), 1)
    print(f'\n{lab}: thin share beyond 1,500 km of subduction {base:.3f}; all {a[sel & thin].sum()/a[sel].sum():.3f}')
    print('  dsub km   | cont area (of sel) | thin share | excess thin area (of total excess) | cum')
    B = list(range(0, 1600, 100)) + [20000]
    ex = []
    for lo, hi in zip(B[:-1], B[1:]):
        k = sel & (dsub >= lo) & (dsub < hi)
        ex.append(max(a[k & thin].sum()-base*a[k].sum(), 0.))
    E = sum(ex); cum = 0.
    for (lo, hi), e in zip(zip(B[:-1], B[1:]), ex):
        k = sel & (dsub >= lo) & (dsub < hi); cum += e
        print(f'  {lo:5d}-{hi:5d} | {a[k].sum()/a[sel].sum():.3f} | {a[k & thin].sum()/max(a[k].sum(),1):.3f} | {e/max(E,1):.3f} | {cum/max(E,1):.2f}')
    w = np.array(ex); mids = np.array([(lo+min(hi, 1600))/2 for lo, hi in zip(B[:-1], B[1:])])
    mu = (w*mids).sum()/w.sum(); sd = np.sqrt((w*(mids-mu)**2).sum()/w.sum())
    print(f'  excess-weighted distance: mean {mu:.0f} km, sd {sd:.0f} km')
