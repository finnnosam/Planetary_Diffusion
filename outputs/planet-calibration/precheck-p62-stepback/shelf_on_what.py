# Pass 62: what the shelf stands on. Coast = land (z > 0) next to sea. Earth: ETOPO on n = 217, CRUST1.0;
# planet: control (combo-01 of sweep-20261101-115244) and combo-02 of sweep-20261101-191306 at 0 Myr, six seeds.
import numpy as np, glob
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, sediment_bulk_m
from planet_tectonics.plates import PARAMS
from planet_earth.reference import ice_mask
g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))
def m(T, S, z, keep):
    cont = T >= 20e3; land = z > 0; sea = ~land
    coast = land & np.any(np.where(ok, sea[np.where(ok, nb, 0)], False), 1) & keep
    dl = gl.distance_to(land, nb, g.dist, 3000e3)/1e3
    shelf = sea & (z > -200) & keep; slope = sea & (z <= -200) & (z > -2000) & keep
    dc = cont & sea & keep
    L = coast.sum()*sp/1e6
    sed_os = S*a*(~cont & sea & (dl <= 150) & keep)
    return [L, a[cont & keep].sum()/1e12, a[dc].sum()/1e12, a[dc].sum()/a[cont & keep].sum(), a[dc].sum()/1e6/(L*1e3),
            np.percentile(dl[dc], 50), np.percentile(dl[dc], 90),
            a[shelf].sum()/1e12, a[shelf & cont].sum()/max(a[shelf].sum(), 1), a[slope].sum()/1e12, a[slope & cont].sum()/max(a[slope].sum(), 1),
            (S*a*dc).sum()/1e15, sed_os.sum()/1e15]
H = 'coast k km | cont | drowned cont | share | per km coast | dist from land p50 p90 | shelf (0..-200) Mkm2 on-cont | slope (-200..-2000) on-cont | sed on drowned cont / ocean side <=150 km (Mkm3)'
print(H)
E = np.load('/home/claude/pre/earth_nodes.npz'); zE = np.load('outputs/planet-reference/earth-n217.npz')['elevation_m'].astype(float)
ice = ice_mask(g)
print('Earth (ice out)', ' '.join('%7.2f' % x for x in m(E['T'], E['S'], zE, ~ice)))
print('Earth (all)    ', ' '.join('%7.2f' % x for x in m(E['T'], E['S'], zE, np.ones(g.size, bool))))
for name, pat in (('control', 'sweep-20261101-115244/combo-01'), ('combo-02', 'sweep-20261101-191306/combo-02')):
    R = []
    for d in sorted(glob.glob(f'/home/claude/pre/sw/{pat}/seed*')):
        z = np.load(f'{d}/snapshots/t+0000.0.npz')
        R.append(m(column_thickness_m(z, PARAMS), sediment_bulk_m(z), z['elevation_m'].astype(float), np.ones(g.size, bool)))
    print('%-15s' % name, ' '.join('%7.2f' % x for x in np.mean(R, 0)))
