"""Pass 96 probe: what the drowned continent is. Label-free continent (crust + bulk sediment >= 20 km) by elevation
about its sea: -1000 ... -200 / -200 ... 0 / 0 ... 500 m, Earth (ETOPO, CRUST1.0, ice sheets out) against the planet
at 0 Myr. Per band: share of the surface; area medians of crystalline crust, bulk sediment, distance from the ocean
(km; distance_to_ocean_km of the label-free continent's complement); the share of the band within 300 km of the
continent's edge; planet only: dyn, lithosphere age, crust age, the share with crust age < 100 Myr, the share whose
crust was thinned (crystalline < 30 km beyond 300 km). Read, not a rule.

Usage (repo root): python3 drowned_band.py <arm folder> [...]
"""
import sys
from pathlib import Path

import numpy as np

from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.coast_budget import wmed
from planet_reviews.coast_profile import distance_to_ocean_km
from planet_reviews.crust_review import sediment_bulk_m

BANDS = ((-1000., -200.), (-200., 0.), (0., 500.), (500., 1e9))


def rows(g, Tc, S, z, extra=None):
    a = g.area
    c = (Tc+S) >= 20e3
    d = distance_to_ocean_km(g, c)
    out = []
    for lo, hi in BANDS:
        m = c & (z >= lo) & (z < hi)
        r = {'share': a[m].sum()/a.sum(), 'Tc': wmed(Tc[m], a[m])/1e3, 'Tc p25': np.nan, 'S': wmed(S[m], a[m]),
             'dist': wmed(d[m], a[m]), 'edge300': a[m & (d < 300)].sum()/max(a[m].sum(), 1e-9)}
        o = np.argsort(Tc[m]); cw = np.cumsum(a[m][o])
        if cw.size:
            r['Tc p25'] = Tc[m][o][np.searchsorted(cw, .25*cw[-1])]/1e3
            r['Tc<30'] = a[m & (Tc < 30e3)].sum()/a[m].sum()
            r['Tc 30-40'] = a[m & (Tc >= 30e3) & (Tc < 40e3)].sum()/a[m].sum()
        for k, v in (extra or {}).items():
            r[k] = wmed(v[m], a[m]) if not k.startswith('share') else a[m & v].sum()/max(a[m].sum(), 1e-9)
        out.append(r)
    return out


def show(name, rs):
    print(name)
    ks = [k for k in rs[0] if k != 'share']
    for (lo, hi), r in zip(BANDS, rs):
        print(f'   {lo:6.0f} .. {hi:<6.0f} share {r["share"]:.3f} | ' + ' | '.join(f'{k} {r[k]:.2f}' if abs(r[k]) < 10 else f'{k} {r[k]:.0f}' for k in ks))


def main(argv):
    g = CubeSphere(217)
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    from planet_earth.reference import ice_mask
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    ice = ice_mask(g)
    show('Earth (CRUST1.0 total >= 20 km, ice sheets out)', rows(g, np.where(ice, 0, T-S), np.where(ice, 0, S), z))
    t = 't+0000.0'
    if argv and argv[0].startswith('--t='):
        t, argv = argv[0][4:], argv[1:]
    print('snapshot', t)
    for arm in map(Path, argv):
        acc = []
        for s in sorted(arm.glob('seed*')):
            f = np.load(s/'snapshots'/f'{t}.npz')
            Tc = f['crust_thickness_m'].astype(float)
            Sb = sediment_bulk_m(f)
            ex = {'dyn': f['dynamic_topography_m'].astype(float), 'L age': f['lithosphere_age_myr'].astype(float),
                  'crust age': f['crust_age_myr'].astype(float), 'share young<100': f['crust_age_myr'] < 100,
                  'share flag': f['continental'].astype(bool), 'thick_m': f['tectonic_thickening_m'].astype(float)}
            acc.append(rows(g, Tc, Sb, f['elevation_m'].astype(float), ex))
        show(f'{arm.name} ({len(acc)} seeds, means)', [{k: float(np.nanmean([a[i][k] for a in acc])) for k in acc[0][i]}
                                                       for i in range(len(BANDS))])


if __name__ == '__main__':
    main(sys.argv[1:])
