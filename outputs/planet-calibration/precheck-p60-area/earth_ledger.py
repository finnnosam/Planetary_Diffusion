# Pass 60: Earth's present continental area ledger. (1) Bird 2003 PB2002 steps by class: area rate = divergence x
# length. (2) GSRM v2.1 (Kreemer et al. 2014) areal dilatation exx+eyy integrated over continental crust (CRUST1.0
# total >= 20 km) on n = 217, split by sign and by the nearest Bird boundary class within 300 km.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_landscape import graph_landscape as gl
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_reviews.margin_width import earth_nodes
r = read_dbf('data/tectonicplates-master/PB2002_steps.dbf')
cl, L, v = r['STEPCLASS'], r['STEPLENGTH'], r['VELOCITYDI']
print('(1) Bird PB2002 steps, area rate (Mkm2/Myr) = divergence x length:')
for c in ('CCB', 'CRB', 'CTF', 'SUB', 'OCB'):
    m = cl == c
    print('   %s %5.1fk km  gain %+.3f  loss %+.3f  net %+.3f' % (c, L[m].sum()/1e3, (np.clip(v[m], 0, None)*L[m]).sum()*1e-6,
          (np.clip(v[m], None, 0)*L[m]).sum()*1e-6, (v[m]*L[m]).sum()*1e-6))
g = CubeSphere(217); a = g.area
T, cont_type = earth_nodes(g)
cont = T >= 20e3
d = np.load('outputs/planet-reference/earth-deformation-n217.npz')
dil = d['exx']+d['eyy']                      # 1/yr
rate = dil*a*1e6/1e12                        # Mkm2/Myr  (1/yr * m2 -> m2/yr; *1e6 yr/Myr; /1e12 m2/Mkm2)
mid = xyz_from_lonlat((r['STARTLONG']+r['FINALLONG'])/2, (r['STARTLAT']+r['FINALLAT'])/2)
near = {}
for c in ('SUB', 'CCB', 'CRB', 'CTF', 'OCB', 'OSR', 'OTF'):
    src = np.zeros(g.size, bool); src[fields.nearest_nodes(g, mid[cl == c])] = True
    near[c] = gl.distance_to(src, g.nbr, g.dist, 3000e3)
D = np.stack([near[c] for c in near]); cls = np.array(list(near))[D.argmin(0)]; dmin = D.min(0)
zone = np.where(dmin <= 300e3, cls, 'far')
print('(2) GSRM areal dilatation over continental crust (CRUST1.0 >= 20 km), Mkm2/Myr, by nearest Bird class <= 300 km:')
print('   %-5s %8s %8s %8s %10s' % ('zone', 'gain', 'loss', 'net', 'area Mkm2'))
for z in ('SUB', 'CCB', 'CRB', 'CTF', 'OCB', 'OSR', 'OTF', 'far', 'all'):
    m = cont & ((zone == z) if z != 'all' else True)
    print('   %-5s %+8.3f %+8.3f %+8.3f %10.1f' % (z, rate[m & (rate > 0)].sum(), rate[m & (rate < 0)].sum(), rate[m].sum(), a[m].sum()/1e12))
m = cont & (zone != 'SUB')
print('   all but SUB-adjacent (interseismic elastic strain at locked megathrusts left out): gain %+.3f loss %+.3f net %+.3f'
      % (rate[m & (rate > 0)].sum(), rate[m & (rate < 0)].sum(), rate[m].sum()))
for lo, hi in ((20e3, 30e3), (30e3, 45e3), (45e3, 1e9)):
    m = cont & (T >= lo) & (T < hi)
    print('   crust %2.0f-%2.0f km: gain %+.3f loss %+.3f net %+.3f' % (lo/1e3, min(hi, 99e3)/1e3, rate[m & (rate > 0)].sum(), rate[m & (rate < 0)].sum(), rate[m].sum()))
