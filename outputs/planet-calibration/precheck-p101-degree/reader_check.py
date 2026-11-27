"""Pass 101: the VD2020 reader's convention, every combination of (sine term on m < 0 or m > 0, Condon-Shortley phase,
sign of the sine term), against Holdt 2022 degree 1-30 on CRUST1.0's oceanic nodes (both measured there), on all nodes,
and coefficient by coefficient against Holdt 1-10 by quadrature. Read. Usage (repo root): python3 reader_check.py"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import DATA
from planet_tectonics import mantle_flow as mf
from planet_tectonics import lithosphere as lt
from degree_split import holdt
from earth_degree_corr import band_corr, LMAX
g = CubeSphere(217); a = g.area
h30 = holdt(g, 'degree_1_to_30.grd'); h10 = holdt(g, 'degree_1_to_10.grd')
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
tn = np.array([l.split() for l in (DATA/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[ii, jj]
oc = np.isin(tn, lt.OCEANIC_TYPES)
rows = np.loadtxt(DATA/'mantle'/'vd2020_sph_coeff_high_accuracy_spot.dat', comments='#')
y = {(int(l), int(m)): v for l, m, v, _ in rows}
A10 = mf.point_alm(g.xyz, h10*a/g.radius**2, LMAX)
def alm(sin_neg, cs, sg):
    A = np.zeros((LMAX+1)*(LMAX+2)//2, complex)
    for (l, m), v in y.items():
        if m < 0 or l > LMAX: continue
        if m == 0: A[mf.alm_index(l, 0, LMAX)] = v; continue
        c, s = (v, y[(l, -m)]) if sin_neg else (y[(l, -m)], v)
        A[mf.alm_index(l, m, LMAX)] = ((-1)**m if cs else 1)*(c-sg*1j*s)/np.sqrt(2)
    return A*1e3
def wc(x, y, w):
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    return np.average((x-xm)*(y-ym), weights=w)/np.sqrt(np.average((x-xm)**2, weights=w)*np.average((y-ym)**2, weights=w))
print(f'oceanic nodes {a[oc].sum()/a.sum():.3f} of the surface')
for sin_neg in (True, False):
    for cs in (False, True):
        for sg in (1, -1):
            V = alm(sin_neg, cs, sg); f = mf.synthesis(V, LMAX, g.xyz)
            print(f'sine on m<0 {sin_neg!s:5} CS {cs!s:5} sine sign {sg:+d}: r(1-30) ocean {wc(f[oc], h30[oc], a[oc]):+.3f} all {wc(f, h30, a):+.3f} | '
                  f'coefficients vs Holdt 1-10: l 1-3 {band_corr(V, A10, LMAX, 1, 3):+.2f}, 4-10 {band_corr(V, A10, LMAX, 4, 10):+.2f}')
print('north-south mirror (x (-1)^(l+m)) of each convention:')
L, M = mf.alm_l(LMAX), mf.alm_m(LMAX)
for sin_neg in (True, False):
    for cs in (False, True):
        for sg in (1, -1):
            V = alm(sin_neg, cs, sg)*(-1.)**(L+M); f = mf.synthesis(V, LMAX, g.xyz)
            print(f'   sine on m<0 {sin_neg!s:5} CS {cs!s:5} sine sign {sg:+d}: r(1-30) ocean {wc(f[oc], h30[oc], a[oc]):+.3f}')
