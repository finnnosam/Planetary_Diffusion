"""Pass 96 probe: where the continental crust sits on the dynamic topography, Earth against the planet.

Earth: CRUST1.0 total thickness (crsthk, sediment included) >= 20 km (the label-free continent), Holdt et al. 2022
degree 1-30 (water-loaded as published; the continental convention unconfirmed, journal 73), ice sheets kept in
the crust mask (the field does not depend on ice). Planet: crust + bulk sediment >= 20 km at 0 Myr, the
air-loaded field (dynamic_topography_air_m: the field before the selection by the sea's water load) and the
applied one. Splits the land's lift into placement (the crust's mean on the field) and selection (land against
drowned crust). Area means, m. Read, not a rule.

Usage (repo root): python3 dyn_crust.py <arm folder> [...]
"""
import sys
from pathlib import Path

import numpy as np

from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import sediment_bulk_m

sys.path.insert(0, str(Path(__file__).parent))
from dyn_inland import holdt   # noqa: E402

BINS = (20e3, 30e3, 35e3, 40e3, 45e3, 55e3, 1e9)


def am(v, w, m):
    return float((v*w)[m].sum()/max(w[m].sum(), 1e-9))


def read(g, T, z, dyn, name):
    a = g.area
    c = T >= 20e3
    land = c & (z > 0)
    wet = c & (z <= 0)
    pos = dyn > 0
    r = {'crust': am(dyn, a, c), 'land': am(dyn, a, land), 'drowned': am(dyn, a, wet), 'ocean': am(dyn, a, ~c),
         'crust on dyn>0': float(a[c & pos].sum()/a[c].sum()),
         'crust share where dyn>0 / <0': (float(a[c & pos].sum()/a[pos].sum()), float(a[c & ~pos].sum()/a[~pos].sum())),
         'bins': [am(dyn, a, c & (T >= lo) & (T < hi)) for lo, hi in zip(BINS[:-1], BINS[1:])]}
    return r


def show(name, r):
    print(f"{name}: mean dyn on the continent {r['crust']:+.0f} m (land {r['land']:+.0f}, drowned {r['drowned']:+.0f}), "
          f"ocean {r['ocean']:+.0f}; continent on dyn > 0 {r['crust on dyn>0']:.2f}; continent's share of the surface "
          f"where dyn > 0 / < 0 {r['crust share where dyn>0 / <0'][0]:.3f} / {r['crust share where dyn>0 / <0'][1]:.3f}")
    print('   by total crust (km) 20-30 / 30-35 / 35-40 / 40-45 / 45-55 / > 55: '
          + ' / '.join(f'{v:+.0f}' for v in r['bins']))


def main(argv):
    g = CubeSphere(217)
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    show('Earth (CRUST1.0 >= 20 km; Holdt 1-30, water-loaded)', read(g, T, z, holdt(g), 'Earth'))
    for arm in map(Path, argv):
        for key in ('dynamic_topography_air_m', 'dynamic_topography_m'):
            rs = []
            for s in sorted(arm.glob('seed*')):
                f = np.load(s/'snapshots'/'t+0000.0.npz')
                Tp = f['crust_thickness_m'].astype(float)+sediment_bulk_m(f)
                rs.append(read(g, Tp, f['elevation_m'].astype(float), f[key].astype(float), arm.name))
            m = {k: np.mean([r[k] for r in rs], axis=0) for k in rs[0]}
            show(f'{arm.name} {key} ({len(rs)} seeds, means)', m)
            print('   continent mean by seed: ' + ' / '.join(f"{r['crust']:+.0f}" for r in rs)
                  + ' | land by seed: ' + ' / '.join(f"{r['land']:+.0f}" for r in rs))


if __name__ == '__main__':
    main(sys.argv[1:])
