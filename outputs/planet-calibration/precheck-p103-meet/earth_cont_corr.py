"""Pass 103 probe: how Earth's measured residual topography meets Earth's continents, by degree, and whether the oceans
constrain it. Earth's field: VD2020's posterior mean (the reader of pass 101, earth_degree_corr.py: cosine on m < 0,
sine on m > 0, Condon-Shortley) and Holdt et al. (2022) degree 1-30 by quadrature (an independent compilation). Earth's
continents: CRUST1.0 on the n = 217 nodes, the label-free continent (crust + sediment >= 20 km) and its normal part
(crust >= 30 km), as indicator functions expanded to l 1-30. Prints:
  1. the degree correlation r_l of each field with each continent function, the 95% band for unrelated fields;
  2. the oceans-only check: per band, the area-weighted correlation over oceanic nodes (where both fields are measured)
     between the field's band and the continent function's band (does the field's long-wavelength part covary with
     the continents where it is constrained);
  3. the field's area mean on the continent by distance from the ocean, per band (the two compilations).
Read, not a rule. Usage (repo root): python3 earth_cont_corr.py"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_reviews.coast_profile import distance_to_ocean_km
P101 = Path('outputs/planet-calibration/precheck-p101-degree').resolve()
sys.path.insert(0, str(P101))
from degree_split import holdt, corr
from earth_degree_corr import vd2020_alm, deg_corr, band_corr, LMAX

g = CubeSphere(217); a = g.area
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
S = np.loadtxt(DATA/'crust1'/'sedthk.xyz')[:, 2].reshape(180, 360)[ii, jj]*1e3
c20 = (T+S) >= 20e3
c30 = c20 & (T >= 30e3)
ocean = ~c20
V = vd2020_alm(False, True)
h30 = holdt(g, 'degree_1_to_30.grd')
H = mf.point_alm(g.xyz, h30*a/g.radius**2, LMAX)
L = mf.alm_l(LMAX)
C = {'c20 (crust+sed >= 20 km)': mf.point_alm(g.xyz, (c20-c20.mean()).astype(float)*a/g.radius**2, LMAX),
     'c30 (crust >= 30 km)': mf.point_alm(g.xyz, (c30-c30.mean()).astype(float)*a/g.radius**2, LMAX)}
F = {'VD2020': V, 'Holdt 1-30': H}
print(f'continent shares: c20 {a[c20].sum()/a.sum():.3f}, c30 {a[c30].sum()/a.sum():.3f}')
print('\n1. degree correlation r_l (95% band for unrelated fields +-1.96/sqrt(2l+1))')
R = {}
for fn, Fa in F.items():
    for cn, Ca in C.items():
        R[(fn, cn)] = deg_corr(Fa, Ca, LMAX)
cols = list(R)
print('   l  | ' + ' | '.join(f'{f[:6]}/{c[:3]}' for f, c in cols) + ' | band')
for l in range(1, 31):
    print(f'   {l:2d} | ' + ' | '.join(f'{R[k][l]:+10.2f}' for k in cols) + f' | +-{1.96/np.sqrt(2*l+1):.2f}')
for lo, hi in ((1, 3), (4, 10), (11, 30), (1, 30)):
    print(f'   band l {lo}-{hi}: ' + ', '.join(f'{f}/{c[:3]} {band_corr(F[f], C[c], LMAX, lo, hi):+.2f}' for f, c in cols))
print(f'   VD2020 vs Holdt by band: ' + ', '.join(f'l {lo}-{hi} {band_corr(V, H, LMAX, lo, hi):+.2f}' for lo, hi in ((1, 3), (4, 10), (11, 30))))
print('\n2. oceans only (oceanic nodes, ~0.6 of the surface): corr(field band, continent-function band)')
for lo, hi in ((1, 3), (4, 10), (11, 30), (1, 30)):
    s = (L >= lo) & (L <= hi)
    out = []
    for fn, Fa in F.items():
        f = mf.synthesis(np.where(s, Fa, 0), LMAX, g.xyz)
        for cn, Ca in C.items():
            c = mf.synthesis(np.where(s, Ca, 0), LMAX, g.xyz)
            out.append(f'{fn}/{cn[:3]} ocean {corr(f[ocean], c[ocean], a[ocean]):+.2f} all {corr(f, c, a):+.2f}')
    print(f'   l {lo}-{hi}: ' + '; '.join(out))
print('\n3. field area mean (m, air-loaded VD2020 x 2300/3300; Holdt as published) on the continent by distance from the ocean')
d = distance_to_ocean_km(g, c20)
bands = ((0, 300), (300, 600), (600, 1200), (1200, 1e9))
for fn, Fa in F.items():
    k = 2300./3300. if fn == 'VD2020' else 1.
    for lo, hi in ((1, 3), (4, 10), (11, 30), (1, 30)):
        f = mf.synthesis(np.where((L >= lo) & (L <= hi), Fa, 0), LMAX, g.xyz)*k
        row = [f'{(a[m]*f[m]).sum()/a[m].sum():+6.0f}' for m in [c20 & (d >= b0) & (d < b1) for b0, b1 in bands]]
        print(f'   {fn:10} l {lo:2d}-{hi:2d}: ocean {(a[ocean]*f[ocean]).sum()/a[ocean].sum():+5.0f} | continent c30 '
              f'{(a[c30]*f[c30]).sum()/a[c30].sum():+5.0f} | by km from the ocean 0-300 / 300-600 / 600-1200 / >=1200: ' + ' / '.join(row))
np.savez(Path(__file__).with_name('earth_cont_corr.npz'), **{f'{f}|{c}': R[(f, c)] for f, c in cols})
