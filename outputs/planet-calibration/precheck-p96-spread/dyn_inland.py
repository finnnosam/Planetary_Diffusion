"""Pass 96 probe: dynamic topography on land by distance inland and by elevation, Earth (Holdt et al. 2022, degree
1-30, as level_split reads it) against the planet's field (dynamic_topography_m), on the same bands as coast_budget.

Question: the planet's interiors stand ~500 m above Earth's against their Airy line (coast_budget, 600+ km), and
~260-390 m of it is the planet's dynamic topography. Does Earth's land show the same lift of its interiors by its
own dynamic topography (the selection 'continents emerge where the field lifts them' acts on Earth too), or is the
planet's placement the defect? Output only: area means and medians by band; read, not a rule.

Usage (repo root): python3 dyn_inland.py <arm folder> [<arm folder> ...]
"""
import sys
from pathlib import Path

import numpy as np
import h5py

from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.coast_profile import BINS_KM, distance_to_ocean_km
from planet_reviews.coast_budget import wmed

ZB = (0., 200., 500., 1000., 2000., 1e9)


def holdt(g):
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    D = DATA/'residual_topography'/'Residual depthtopographic measurements'/'Residual depthtopographic measurements'/'Holdt_et_al_2022_SupportingInformation'
    with h5py.File(D/'degree_1_to_30.grd', 'r') as f:
        lo, la, v = f['lon'][:], f['lat'][:], f['z'][:]
    i = np.clip(np.round((lat-la[0])/(la[1]-la[0])).astype(int), 0, la.size-1)
    j = np.round(((lon % 360)-lo[0])/(lo[1]-lo[0])).astype(int) % lo.size
    return np.nan_to_num(v[i, j]*1e3)


def table(g, z, dyn, land):
    a = g.area
    d = distance_to_ocean_km(g, land)
    out = {'all land': float((a*dyn)[land].sum()/a[land].sum()),
           'ocean': float((a*dyn)[~land].sum()/a[~land].sum())}
    dist = []
    for lo, hi in zip(BINS_KM[:-1], BINS_KM[1:]):
        m = land & (d >= lo) & (d < hi)
        dist.append((a[m].sum()/a[land].sum(), float((a*dyn)[m].sum()/max(a[m].sum(), 1)), wmed(dyn[m], a[m])))
    elev = []
    for lo, hi in zip(ZB[:-1], ZB[1:]):
        m = land & (z >= lo) & (z < hi)
        elev.append((a[m].sum()/a.sum(), float((a*dyn)[m].sum()/max(a[m].sum(), 1))))
    # the share of land on positive dyn, and land's share of the surface where dyn > 0 / < 0
    pos = dyn > 0
    out['land on dyn>0'] = float(a[land & pos].sum()/a[land].sum())
    out['land share where dyn>0 / <0'] = (float(a[land & pos].sum()/a[pos].sum()), float(a[land & ~pos].sum()/a[~pos].sum()))
    out['corr(z, dyn) on land'] = float(np.corrcoef(z[land], dyn[land])[0, 1])
    return out, dist, elev


def show(name, res):
    out, dist, elev = res
    print(f'{name}: mean dyn on land {out["all land"]:+.0f} m, on sea {out["ocean"]:+.0f}; land on dyn > 0 {out["land on dyn>0"]:.2f}; '
          f'land share where dyn > 0 / < 0 {out["land share where dyn>0 / <0"][0]:.3f} / {out["land share where dyn>0 / <0"][1]:.3f}; '
          f'corr(z, dyn) on land {out["corr(z, dyn) on land"]:.2f}')
    print('   by distance inland (km): land share | dyn mean | dyn median')
    for (lo, hi), r in zip(zip(BINS_KM[:-1], BINS_KM[1:]), dist):
        print(f'     {lo:5.0f}-{hi:<6.0f} {r[0]:.3f} | {r[1]:+5.0f} | {r[2]:+5.0f}')
    print('   by elevation above the sea (m): surface share | dyn mean')
    for (lo, hi), r in zip(zip(ZB[:-1], ZB[1:]), elev):
        print(f'     {lo:5.0f}-{hi:<6.0f} {r[0]:.3f} | {r[1]:+5.0f}')
    print(f'   coast (0-50) to interior (600-1200): dyn median {dist[4][2]-dist[0][2]:+.0f} m, mean {dist[4][1]-dist[0][1]:+.0f} m')


def main(argv):
    g = CubeSphere(217)
    from planet_earth.reference import ice_mask
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    land = (z > 0) & ~ice_mask(g)
    show('Earth (ETOPO, ice sheets out; Holdt 1-30)', table(g, z, holdt(g), land))
    for arm in map(Path, argv):
        acc = []
        for s in sorted(arm.glob('seed*')):
            f = np.load(s/'snapshots'/'t+0000.0.npz')
            e = f['elevation_m'].astype(float)
            acc.append(table(g, e, f['dynamic_topography_m'].astype(float), e > 0))
        o = {k: (np.mean([a[0][k] for a in acc], axis=0)) for k in acc[0][0]}
        o['land share where dyn>0 / <0'] = tuple(o['land share where dyn>0 / <0'])
        dist = np.nanmean([a[1] for a in acc], axis=0); elev = np.nanmean([a[2] for a in acc], axis=0)
        show(f'{arm.name} ({len(acc)} seeds, 0 Myr, means)', (o, dist, elev))
        c = [a[1][4][2]-a[1][0][2] for a in acc]
        print('   coast to interior dyn median by seed: ' + ' / '.join(f'{v:+.0f}' for v in c))


if __name__ == '__main__':
    main(sys.argv[1:])
