"""Pass 96 probe: does the stable cover lie where its column has room for it? Earth against the planet's start.

Domain (Earth): the start cover's own (tests/test_start_sediment.py): CRUST1.0 total >= 20 km, Precambrian /
platform types, > 500 km from Bird's boundaries. The planet: the start's continent (crust + bulk sediment >= 20 km)
at -250 Myr, and the same at 0 Myr. On each: the normal scores (area-weighted) of the bulk cover S, the crystalline
crust Tc and the bare level b = z - S (1 - rho_sediment / rho_mantle) (the surface on the datum less the cover's own
isostatic lift, the planet's Airy densities; on Earth b holds its C, dyn and residual). Correlations, the partial
correlation of S and b given Tc, and the share of the domain under the datum with its mean cover. Read, not a rule.

Usage (repo root): python3 cover_accommodation.py <arm folder> [...]
"""
import sys
from pathlib import Path

import numpy as np

from planet_grid.paths import OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import sediment_bulk_m
from planet_tectonics.plates import PARAMS
from planet_tectonics.start import _nscore


def corr(a, x, y):
    return float(np.sum(a*x*y)/np.sqrt(np.sum(a*x*x)*np.sum(a*y*y)))


def read(a, S, Tc, z, p):
    b = z-S*(1-p['rho_sediment']/p['rho_mantle'])
    zs, zt, zb, zz = _nscore(S, a), _nscore(Tc, a), _nscore(b, a), _nscore(z, a)
    r_sb, r_st, r_tb = corr(a, zs, zb), corr(a, zs, zt), corr(a, zt, zb)
    part = (r_sb-r_st*r_tb)/np.sqrt((1-r_st**2)*(1-r_tb**2))
    under = z < 0
    return {'r(S,b)': r_sb, 'r(S,Tc)': r_st, 'r(Tc,b)': r_tb, 'r(S,b|Tc)': float(part), 'r(S,z)': corr(a, zs, zz),
            'under': float(a[under].sum()/a.sum()), 'S under': float((a*S)[under].sum()/max(a[under].sum(), 1e-9)),
            'S above': float((a*S)[~under].sum()/a[~under].sum()), 'S mean': float((a*S).sum()/a.sum())}


def show(name, r):
    print(f"{name}: r(S, b) {r['r(S,b)']:+.2f}, r(S, Tc) {r['r(S,Tc)']:+.2f}, r(Tc, b) {r['r(Tc,b)']:+.2f}, "
          f"partial r(S, b | Tc) {r['r(S,b|Tc)']:+.2f}, r(S, z) {r['r(S,z)']:+.2f}; under the datum {r['under']:.3f} "
          f"(cover there {r['S under']:.0f} m, above {r['S above']:.0f}, mean {r['S mean']:.0f} m bulk)")


def main(argv):
    g = CubeSphere(217); a = g.area
    from planet_reviews.interior_probe import GROUPS
    from planet_earth.bird import earth_boundaries
    from planet_earth.crust1 import C1, crust1_cell
    from planet_earth.reference import ice_mask
    cell = crust1_cell(g)
    T = np.loadtxt(C1/'crsthk.xyz')[:, 2][cell]*1e3; S = np.loadtxt(C1/'sedthk.xyz')[:, 2][cell]*1e3
    ty = np.array([l.split() for l in (C1/'CNtype1-1.txt').read_text().splitlines() if l.strip()]).ravel()[cell]
    dbnd = gl.distance_to(earth_boundaries(g), g.nbr, g.dist, 5000e3)/1e3
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    m = (T >= 20e3) & np.isin(ty, GROUPS['Precambrian/platform']) & (dbnd > 500) & ~ice_mask(g)
    show('Earth stable domain (ice out)', read(a[m], S[m], (T-S)[m], z[m], PARAMS))
    m2 = (T >= 20e3) & ~ice_mask(g)
    show('Earth, all continent >= 20 km (ice out)', read(a[m2], S[m2], (T-S)[m2], z[m2], PARAMS))
    for arm in map(Path, argv):
        for t in ('t-0250.0', 't+0000.0'):
            rs = []
            for s in sorted(arm.glob('seed*')):
                f = np.load(s/'snapshots'/f'{t}.npz')
                Tc = f['crust_thickness_m'].astype(float); Sb = sediment_bulk_m(f)
                c = (Tc+Sb) >= 20e3
                zd = f['elevation_m'].astype(float)+float(f['sea_level_m'])
                rs.append(read(a[c], Sb[c], Tc[c], zd[c], PARAMS))
            show(f'{arm.name} {t} ({len(rs)} seeds, means)', {k: float(np.mean([r[k] for r in rs])) for k in rs[0]})


if __name__ == '__main__':
    main(sys.argv[1:])
