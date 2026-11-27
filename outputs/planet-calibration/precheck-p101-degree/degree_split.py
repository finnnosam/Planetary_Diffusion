"""Pass 101 probe (journal 136): the placed dynamic field split by spherical-harmonic degree on the sets that carry the
target's two failures: the drowned normal crust (label-free continent, crystalline >= 30 km, -1000 ... -200 m about
its sea), the drowned thin crust, and land by distance from the ocean (0-300 / 300-600 / 600-1200 / >= 1200 km).
Planet: the whole air-loaded field (dynamic_topography_air_m, swells included) by quadrature, synthesised by band
l 1-3 / 4-10 / 11-31 and water-loaded under the sea as the model does (x rho_m / (rho_m - rho_w)); 'placed' is the
flow's coefficients alone (flow_topography_alm_m; with flow_placed 'earth_corr' only the placed share), 'rest' the
remainder of dynamic_topography_m (swells, and under 'earth_corr' the plate-blind part). Earth: Holdt et al. (2022) degree 1-10 and 1-30 grids (11-30
= their difference), ETOPO and CRUST1.0 on the same nodes, ice sheets out. Also per band: the area-weighted correlation
of the field with the distance to young slabs (the design check P-b by degree) and the interior-minus-coast contrast.
Read, not a rule. Usage (repo root): python3 degree_split.py <arm> [--tags t+0000.0,t-0020.0,t-0040.0]"""
import sys
from pathlib import Path
import numpy as np
import h5py
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.coast_profile import distance_to_ocean_km
from planet_reviews.crust_review import sediment_bulk_m, run_params
from planet_tectonics import mantle_flow as mf

SEEDS = range(1001, 1007)
SETS = ('continent Tc>=30', 'drowned Tc>=30', 'drowned Tc<30', 'land 0-300', 'land 300-600', 'land 600-1200',
        'land >=1200', 'all land')
HOLDT = DATA/'residual_topography'/'Residual depthtopographic measurements'/'Residual depthtopographic measurements'/'Holdt_et_al_2022_SupportingInformation'


def masks(g, Tc, Sb, zs):
    c = (Tc+Sb) >= 20e3
    land = c & (zs >= 0)
    d = distance_to_ocean_km(g, land)
    dr = c & (zs >= -1000) & (zs < -200)
    return [c & (Tc >= 30e3), dr & (Tc >= 30e3), dr & (Tc < 30e3), land & (d < 300), land & (d >= 300) & (d < 600),
            land & (d >= 600) & (d < 1200), land & (d >= 1200), land], (land & (d < 50), land & (d >= 600) & (d < 1200))


def means(a, fields, ms):
    return np.array([[float((a[m]*f[m]).sum()/max(a[m].sum(), 1e-30)) for f in fields] for m in ms])


def corr(x, y, w):
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.average((x-xm)*(y-ym), weights=w)/np.sqrt(np.average((x-xm)**2, weights=w)*np.average((y-ym)**2, weights=w)))


def dist_to(g, src):
    s = g.xyz[src]; d = np.empty(g.size)
    for i in range(0, g.size, 4000):
        d[i:i+4000] = np.arccos(np.clip((g.xyz[i:i+4000]@s.T).max(1), -1, 1))*g.radius/1e3
    return d


def holdt(g, name):
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    with h5py.File(HOLDT/name, 'r') as f:
        lo, la, v = f['lon'][:], f['lat'][:], f['z'][:]
    i = np.clip(np.round((lat-la[0])/(la[1]-la[0])).astype(int), 0, la.size-1)
    j = np.round(((lon % 360)-lo[0])/(lo[1]-lo[0])).astype(int) % lo.size
    return np.nan_to_num(v[i, j]*1e3)


def show(title, cols, M, extra=''):
    print(title + extra)
    print('   set               | share  | ' + ' | '.join(f'{c:>7}' for c in cols))
    for s, r in zip(SETS, M):
        print(f'   {s:17} | {r[0]:.4f} | ' + ' | '.join(f'{v:+7.0f}' for v in r[1:]))


def main(argv):
    tags = ['t+0000.0']
    if '--tags' in argv:
        i = argv.index('--tags'); tags = argv[i+1].split(','); argv = argv[:i]+argv[i+2:]
    g = CubeSphere(217); a = g.area
    # Earth
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    from planet_earth.reference import ice_mask
    ice = ice_mask(g)
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    h10, h30 = holdt(g, 'degree_1_to_10.grd'), holdt(g, 'degree_1_to_30.grd')
    ms, (coast, inner) = masks(g, np.where(ice, 0, T-S), np.where(ice, 0, S), z)
    shares = np.array([a[m].sum()/a.sum() for m in ms])
    M = means(a, [h10, h30-h10, h30], ms)
    show('Earth (Holdt 2022; CRUST1.0, ETOPO, ice out)', ['1-10', '11-30', '1-30'], np.c_[shares, M])
    ci = [float((a[inner]*f[inner]).sum()/a[inner].sum()-(a[coast]*f[coast]).sum()/a[coast].sum()) for f in (h10, h30-h10, h30)]
    print('   interior 600-1200 minus coast 0-50: ' + ' / '.join(f'{v:+.0f}' for v in ci) + ' m')
    # planet
    for arm in map(Path, argv):
        R, CI, CO, SH = [], [], [], []
        for s in SEEDS:
            p = run_params(arm/f'seed{s}'); fw_wet = p['rho_mantle']/(p['rho_mantle']-p['rho_water'])
            for t in tags:
                q = np.load(arm/f'seed{s}'/'snapshots'/f'{t}.npz')
                H = q['flow_topography_alm_m'].astype(np.complex128); lm = int(q['flow_lmax']); L = mf.alm_l(lm)
                zs = q['elevation_m'].astype(float)
                fw = np.where(zs < 0, fw_wet, 1.)
                At = mf.point_alm(g.xyz, q['dynamic_topography_air_m'].astype(float)*a/g.radius**2, lm)
                bands = [mf.synthesis(np.where((L >= lo) & (L <= hi), At, 0), lm, g.xyz) for lo, hi in ((1, 3), (4, 10), (11, lm))]
                tot = mf.synthesis(H, lm, g.xyz)
                dyn = q['dynamic_topography_m'].astype(float)
                fl = [b*fw for b in bands]+[tot*fw, dyn-tot*fw, dyn]
                m2, (cst, inn) = masks(g, q['crust_thickness_m'].astype(float), sediment_bulk_m(q), zs)
                R.append(np.c_[[a[m].sum()/a.sum() for m in m2], means(a, fl, m2)])
                CI.append([float((a[inn]*f[inn]).sum()/a[inn].sum()-(a[cst]*f[cst]).sum()/a[cst].sum()) for f in fl])
                young = q['slab_mass_young_kg'] > 0
                dy = dist_to(g, np.flatnonzero(young)) if young.any() else None
                CO.append([corr(b, dy, a) if dy is not None else np.nan for b in bands+[sum(bands)]])
                SH.append(float(np.sqrt((a*(dyn-tot*fw)**2).sum()/a.sum())))
        show(f'{arm.name} ({len(SEEDS)} seeds x {len(tags)} snapshots {",".join(tags)}; means)',
             ['1-3', '4-10', '11-31', 'placed', 'rest', 'dyn'], np.mean(R, 0))
        print('   interior 600-1200 minus coast 0-50: ' + ' / '.join(f'{v:+.0f}' for v in np.mean(CI, 0)) + ' m')
        print('   corr(field, distance to young slabs) by band 1-3 / 4-10 / 11-31 / all: ' +
              ' / '.join(f'{v:+.2f}' for v in np.nanmean(CO, 0)) + f'  (by seed, all: {", ".join(f"{c[-1]:+.2f}" for c in CO)})')
        print(f'   rms of the rest (dyn minus the placed flow): {np.mean(SH):.0f} m')


if __name__ == '__main__':
    main(sys.argv[1:])
