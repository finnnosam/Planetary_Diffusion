import numpy as np, sys
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf, plates as PL
from planet_reviews.crust_review import sediment_bulk_m
from planet_reviews.coast_profile import distance_to_ocean_km
g = CubeSphere(217); a = g.area
def dist_to(src):
    s = g.xyz[src]; d = np.empty(g.size)
    for i in range(0, g.size, 4000):
        d[i:i+4000] = np.arccos(np.clip((g.xyz[i:i+4000]@s.T).max(1), -1, 1))*g.radius/1e3
    return d
def corr(x, y):
    xm, ym = np.average(x, weights=a), np.average(y, weights=a)
    return np.average((x-xm)*(y-ym), weights=a)/np.sqrt(np.average((x-xm)**2, weights=a)*np.average((y-ym)**2, weights=a))
for arm in sys.argv[1:]:
    q = np.load(f'{arm}/seed1026/snapshots/t+0000.0.npz')
    air = q['dynamic_topography_air_m'].astype(float); zs = q['elevation_m'].astype(float)
    fw = np.where(zs < 0, 3300/(3300-1030), 1.); L = mf.alm_l(31)
    A = mf.point_alm(g.xyz, air*a/g.radius**2, 31)
    B = [mf.synthesis(np.where((L >= lo) & (L <= hi), A, 0), 31, g.xyz) for lo, hi in ((1, 3), (4, 10), (11, 31))]
    Tc = q['crust_thickness_m'].astype(float); Sb = sediment_bulk_m(q); c = (Tc+Sb) >= 20e3
    land = c & (zs >= 0); d = distance_to_ocean_km(g, land); ds = dist_to(np.flatnonzero(q['boundary_type'] == PL.SUBDUCTION))
    sets = {'drowned Tc>=30': c & (Tc >= 30e3) & (zs >= -1000) & (zs < -200), 'land 600-1200': land & (d >= 600) & (d < 1200),
            'land >=1200': land & (d >= 1200), 'cont30 0-600 km from subduction': c & (Tc >= 30e3) & (ds < 600)}
    print(arm.split('/')[-1], '(seed 1026, 0 Myr; whole air field by band, water-loaded under the sea)')
    for k, m in sets.items():
        W = a[m].sum()
        print(f'   {k:32} share {W/a.sum():.4f}: ' + ' / '.join(f'{(a[m]*b[m]*fw[m]).sum()/max(W,1e-30):+.0f}' for b in B) + ' m (l 1-3 / 4-10 / 11-31)')
    print('   corr(band, distance to subduction): ' + ' / '.join(f'{corr(b, ds):+.2f}' for b in B) + f'; whole {corr(air, ds):+.2f}')
