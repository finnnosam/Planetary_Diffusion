"""Pass 101 probe (journal 136): the continent (label-free, crystalline >= 30 km) by distance from subduction, Earth
against the planet. Per band (0-300 / 300-600 / 600-1000 / 1000-2000 / >= 2000 km): share of the surface, the dynamic
field's area mean at degree 1-10 and above 10 (Earth: Holdt 2022 grids; planet: flow coefficients water-loaded under the
sea), the share drowned (-1000 ... -200 m about the sea) and the area median elevation about the sea. Earth's
subduction: the midpoints of Bird's (2003) SUB steps; the planet's: nodes with boundary_type SUBDUCTION. Separates two
readings of the drowned normal crust: more continent near trenches (geography) or deeper lows at the same distance.
Read, not a rule. Usage (repo root): python3 trench_field.py <arm>"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import sediment_bulk_m, run_params
from planet_reviews.trench_origin import read_dbf, _xyz
from planet_tectonics import mantle_flow as mf
from planet_tectonics import plates as PL
sys.path.insert(0, str(Path(__file__).resolve().parent))
from degree_split import holdt, dist_to

BANDS = (0., 300., 600., 1000., 2000., 1e9)


def wmed(x, w):
    if x.size == 0:
        return np.nan
    o = np.argsort(x); c = np.cumsum(w[o]); return float(x[o][np.searchsorted(c, .5*c[-1])])


def rows(g, Tc, Sb, zs, ds, lo10, hi10):
    a = g.area; c = ((Tc+Sb) >= 20e3) & (Tc >= 30e3); out = []
    for b0, b1 in zip(BANDS[:-1], BANDS[1:]):
        m = c & (ds >= b0) & (ds < b1); W = max(a[m].sum(), 1e-30)
        out.append([a[m].sum()/a.sum(), (a[m]*lo10[m]).sum()/W, (a[m]*hi10[m]).sum()/W,
                    a[m & (zs >= -1000) & (zs < -200)].sum()/W, a[m & (zs >= 0)].sum()/W, wmed(zs[m], a[m])])
    return np.array(out)


def show(title, M):
    print(title)
    print('   distance from subduction | share  | dyn 1-10 | dyn >10 | drowned | land  | z p50')
    for (b0, b1), r in zip(zip(BANDS[:-1], BANDS[1:]), M):
        print(f'   {b0:5.0f} - {b1:<6.0f}           | {r[0]:.4f} | {r[1]:+8.0f} | {r[2]:+7.0f} | {r[3]:.3f}   | {r[4]:.3f} | {r[5]:+5.0f}')


def main(argv):
    g = CubeSphere(217)
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    from planet_earth.reference import ice_mask
    ice = ice_mask(g)
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    h10, h30 = holdt(g, 'degree_1_to_10.grd'), holdt(g, 'degree_1_to_30.grd')
    s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
    sub = np.asarray(s['STEPCLASS']) == 'SUB'
    m = _xyz(np.asarray(s['STARTLONG'])[sub], np.asarray(s['STARTLAT'])[sub])+_xyz(np.asarray(s['FINALLONG'])[sub], np.asarray(s['FINALLAT'])[sub])
    m /= np.linalg.norm(m, axis=1, keepdims=True)
    d = np.empty(g.size)
    for i in range(0, g.size, 4000):
        d[i:i+4000] = np.arccos(np.clip((g.xyz[i:i+4000]@m.T).max(1), -1, 1))*g.radius/1e3
    show('Earth (CRUST1.0, ETOPO, ice out; Holdt 1-10 / 11-30; Bird SUB)',
         rows(g, np.where(ice, 0, T-S), np.where(ice, 0, S), z, d, h10, h30-h10))
    for arm in map(Path, argv):
        R = []
        for sd in range(1001, 1007):
            p = run_params(arm/f'seed{sd}'); fw_wet = p['rho_mantle']/(p['rho_mantle']-p['rho_water'])
            q = np.load(arm/f'seed{sd}'/'snapshots'/'t+0000.0.npz')
            H = q['flow_topography_alm_m'].astype(np.complex128); lm = int(q['flow_lmax']); L = mf.alm_l(lm)
            zs = q['elevation_m'].astype(float); fw = np.where(zs < 0, fw_wet, 1.)
            lo = mf.synthesis(np.where(L <= 10, H, 0), lm, g.xyz)*fw; hi = mf.synthesis(np.where(L > 10, H, 0), lm, g.xyz)*fw
            ds = dist_to(g, np.flatnonzero(q['boundary_type'] == PL.SUBDUCTION))
            R.append(rows(g, q['crust_thickness_m'].astype(float), sediment_bulk_m(q), zs, ds, lo, hi))
        show(f'{arm.name} (6 seeds, 0 Myr, means)', np.nanmean(R, 0))


if __name__ == '__main__':
    main(sys.argv[1:])
