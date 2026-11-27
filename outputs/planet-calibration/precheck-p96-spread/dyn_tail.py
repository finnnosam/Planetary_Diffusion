"""Pass 96 probe: how much of the continents' drowned low tail is the dynamic topography.

Continents = crust >= 30 km (continent.py's set; planet: crystalline crust, as continent.py reads it), elevation on
the datum (planet: elevation_m + sea_level_m). Quantiles p5 / p10 / p25 / p50 / p75 / p90 and the share under the
datum of z, and of z - dyn (planet: the applied field dynamic_topography_m; Earth: Holdt 2022 l 1-30), and the share
under the datum by total-crust bin. z - dyn is a counterfactual without re-solving isostasy, flexure or the sea: read
as the field's direct share only. Area-weighted. Read, not a rule.

Usage (repo root): python3 dyn_tail.py <arm folder> [...]
"""
import sys
from pathlib import Path

import numpy as np

from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import sediment_bulk_m

sys.path.insert(0, str(Path(__file__).parent))
from dyn_inland import holdt   # noqa: E402

Q = (5, 10, 25, 50, 75, 90)
BINS = (30e3, 35e3, 40e3, 45e3, 55e3)


def wq(v, w, q):
    o = np.argsort(v); c = np.cumsum(w[o]); c = c/c[-1]
    return [float(v[o][np.searchsorted(c, x/100)]) for x in q]


def read(g, Tc, Ttot, z, dyn):
    a = g.area
    m = Tc >= 30e3
    r = {}
    for k, v in (('z', z), ('z-dyn', z-dyn)):
        r[k] = wq(v[m], a[m], Q)+[float(a[m & (v < 0)].sum()/a[m].sum())]
        r[k+' bins'] = [float(a[mb & (v < 0)].sum()/max(a[mb].sum(), 1e-9))
                        for lo, hi in zip(BINS[:-1], BINS[1:]) for mb in [(Ttot >= lo) & (Ttot < hi)]]
    r['dyn under the datum'] = float((a*dyn)[m & (z < 0)].sum()/a[m & (z < 0)].sum())
    r['dyn above'] = float((a*dyn)[m & (z >= 0)].sum()/a[m & (z >= 0)].sum())
    return r


def show(name, r):
    print(name)
    for k in ('z', 'z-dyn'):
        print(f"   {k:6s}: p5/p10/p25/p50/p75/p90 " + ' / '.join(f'{v:.0f}' for v in r[k][:-1])
              + f"; under the datum {r[k][-1]:.3f}; by total crust 30-35/35-40/40-45/45-55 under the datum "
              + ' / '.join(f'{v:.2f}' for v in r[k+' bins']))
    print(f"   mean dyn on the continent under / above the datum {r['dyn under the datum']:+.0f} / {r['dyn above']:+.0f} m")


def main(argv):
    g = CubeSphere(217)
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    from planet_earth.reference import ice_mask
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    keep = ~ice_mask(g)
    Tc = np.where(keep, T-S, 0.)
    show('Earth (crystalline >= 30 km, ice sheets out; Holdt 1-30)', read(g, Tc, T, z, holdt(g)))
    for arm in map(Path, argv):
        rs = []
        for s in sorted(arm.glob('seed*')):
            f = np.load(s/'snapshots'/'t+0000.0.npz')
            Tc = f['crust_thickness_m'].astype(float)
            zd = f['elevation_m'].astype(float)+float(f['sea_level_m'])
            rs.append(read(g, Tc, Tc+sediment_bulk_m(f), zd, f['dynamic_topography_m'].astype(float)))
        show(f'{arm.name} ({len(rs)} seeds, means)', {k: np.mean([r[k] for r in rs], axis=0) for k in rs[0]})
        print('   under the datum by seed, z / z-dyn: ' + ' / '.join(f"{r['z'][-1]:.2f}-{r['z-dyn'][-1]:.2f}" for r in rs))


if __name__ == '__main__':
    main(sys.argv[1:])
