# Pass 60: ocean-floor sediment (crust + sediment < 20 km) by distance from the label-free continent and from
# trenches; Earth GlobSed (valid ocean nodes; CRUST1.0 sedthk elsewhere) against combo-02 / combo-04 at 0 Myr.
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_landscape import graph_landscape as gl
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_earth.crust1 import C1
from planet_reviews.crust_review import column_thickness_m, run_params, sediment_bulk_m
g = CubeSphere(217); a = g.area; nb = g.nbr
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
S1 = (np.loadtxt(C1/'sedthk.xyz')[:, 2].reshape(180, 360)*1e3)[r, c]
T = np.load('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz')['T'].astype(float)
gs = np.load('outputs/planet-reference/globsed-n217.npz'); G = np.where((gs['valid_share'] > .5) & np.isfinite(gs['sediment_m']), gs['sediment_m'], S1)
b = read_dbf('data/tectonicplates-master/PB2002_steps.dbf')
mid = xyz_from_lonlat((b['STARTLONG']+b['FINALLONG'])/2, (b['STARTLAT']+b['FINALLAT'])/2)
src = np.zeros(g.size, bool); src[fields.nearest_nodes(g, mid[b['STEPCLASS'] == 'SUB'])] = True
D = [(0, 200), (200, 500), (500, 1000), (1000, 1e5)]
def split(sed, cont, trench):
    oc = ~cont; dc = gl.distance_to(cont, nb, g.dist, 1.1e8)/1e3; dt = gl.distance_to(trench, nb, g.dist, 3e6)/1e3
    V = lambda m: (sed[m]*a[m]).sum()/1e15
    return [V(oc)] + [V(oc & (dc >= lo) & (dc < hi)) for lo, hi in D] + [V(oc & (dt < 150)), a[oc & (dt < 150)].sum()/1e12]
hdr = 'ocean total | from continent ' + ' '.join('%d-%d' % (lo, min(hi, 9999)) for lo, hi in D) + ' km | within 150 km of a trench (Mkm3, area Mkm2)'
print(hdr)
print('Earth     ' + ' '.join('%7.0f' % x for x in split(G, T >= 20e3, src)))
for combo in ('combo-02', 'combo-04'):
    R = []
    for seed in range(1001, 1007):
        d = f'outputs/planet-calibration/sweep-20261101-191306/{combo}/seed{seed}'; z = np.load(f'{d}/snapshots/t+0000.0.npz')
        Tp = column_thickness_m(z, run_params(d)); bt = z['boundary_type']
        R.append(split(sediment_bulk_m(z), Tp >= 20e3, bt == 2))
    print('%-9s ' % combo + ' '.join('%7.0f' % x for x in np.mean(R, 0)))
