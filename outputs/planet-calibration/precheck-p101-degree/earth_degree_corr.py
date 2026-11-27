"""Pass 101 probe (journal 136): how much of Earth's measured residual topography the slab flow explains, by degree.
Earth's measured field: Valentine & Davies (2020) posterior-mean coefficients (data/mantle/vd2020_sph_coeff_high_
accuracy_spot.dat, real orthonormal harmonics, km, water-loaded). Earth's slab field: E3 (pass 49, journal 70):
Muller et al. (2019) subduction history through D's unchanged record and kernels, present day, air-loaded
(precheck-p49-e3/e3_H.npz; with today's slab_alpha 'thermo' from e3_thermo.py). The reader's convention (which sign of m is the sine term, Condon-Shortley phase) is
chosen as the one whose synthesis agrees with Holdt et al. (2022) degree 1-30 grid, an independent compilation: a check
of the reader, not a fit. Prints the degree correlation r_l of VD2020 with E3, the 95% band of r_l for two unrelated
fields with 2l+1 coefficients (+-1.96 / sqrt(2l+1)), and the same for the tomographic slab proxy (E1) as a cross-check.
Read, not a rule. Usage (repo root): python3 earth_degree_corr.py"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
sys.path.insert(0, str(Path(__file__).resolve().parent))
from degree_split import holdt

LMAX = 31


def vd2020_alm(sin_neg=True, cs=False):
    rows = np.loadtxt(DATA/'mantle'/'vd2020_sph_coeff_high_accuracy_spot.dat', comments='#')
    a = np.zeros((LMAX+1)*(LMAX+2)//2, np.complex128)
    y = {(int(l), int(m)): v for l, m, v, _ in rows}
    for (l, m), v in y.items():
        if m < 0 or l > LMAX:
            continue
        if m == 0:
            a[mf.alm_index(l, 0, LMAX)] = v
            continue
        c, s = (v, y[(l, -m)]) if sin_neg else (y[(l, -m)], v)
        ph = (-1)**m if cs else 1.
        a[mf.alm_index(l, m, LMAX)] = ph*(c-1j*s)/np.sqrt(2)
    return a*1e3


def deg_corr(A, B, lmax):
    l, m = mf.alm_l(lmax), mf.alm_m(lmax); w = np.where(m == 0, 1., 2.)
    num = np.bincount(l, w*(A*np.conj(B)).real, lmax+1)
    den = np.sqrt(np.bincount(l, w*np.abs(A)**2, lmax+1)*np.bincount(l, w*np.abs(B)**2, lmax+1))
    return np.where(den > 0, num/np.where(den > 0, den, 1), np.nan)


def band_corr(A, B, lmax, lo, hi):
    l, m = mf.alm_l(lmax), mf.alm_m(lmax); w = np.where(m == 0, 1., 2.); s = (l >= lo) & (l <= hi)
    return float((w[s]*(A[s]*np.conj(B[s])).real).sum()/np.sqrt((w[s]*np.abs(A[s])**2).sum()*(w[s]*np.abs(B[s])**2).sum()))


def main():
    g = CubeSphere(217); a = g.area
    h30 = holdt(g, 'degree_1_to_30.grd')
    best = None
    print('reader check: corr(VD2020 synthesis, Holdt 1-30) over all nodes, by convention')
    for sin_neg in (True, False):
        for cs in (False, True):
            f = mf.synthesis(vd2020_alm(sin_neg, cs), LMAX, g.xyz)
            fm, hm = np.average(f, weights=a), np.average(h30, weights=a)
            r = np.average((f-fm)*(h30-hm), weights=a)/np.sqrt(np.average((f-fm)**2, weights=a)*np.average((h30-hm)**2, weights=a))
            print(f'   sine on m<0 {sin_neg!s:5}  Condon-Shortley {cs!s:5}: r {r:+.3f}; rms {np.sqrt(np.average(f**2, weights=a)):.0f} m')
            if best is None or r > best[0]:
                best = (r, sin_neg, cs)
    V = vd2020_alm(best[1], best[2])
    E3 = np.load(Path(__file__).with_name('e3_thermo_H.npz'))           # e3_thermo.py: today's slab_alpha 'thermo'
    H3 = E3['H']; assert int(E3['lmax']) == LMAX
    rows = [('E3 (Earth slab history, D, slab_alpha thermo)', H3),
            ('E3 as run in pass 49 (slab_alpha depth)', np.load(OUTPUTS/'planet-calibration'/'precheck-p49-e3'/'e3_H.npz')['H'])]
    try:
        from planet_reviews import mantle_flow_earth as mfe
        l1, H1, _ = mfe.earth_flow('smean2', True, 300.)
        if l1 == LMAX:
            rows.append(('E1 (SMEAN2 fast > 300 km)', H1))
    except Exception as e:                                    # E1 needs the kernels; optional
        print('   (E1 skipped:', e, ')')
    P = mf.spectrum(V, LMAX)*4*np.pi
    for name, H in rows:
        r = deg_corr(V, H, LMAX)
        print(f'\nVD2020 mean against {name}: degree correlation r_l (95% band for unrelated fields +-1.96/sqrt(2l+1))')
        for l in range(1, 31):
            print(f'   l {l:2d}  r {r[l]:+.2f}  (+-{1.96/np.sqrt(2*l+1):.2f})  VD2020 mean power {P[l]/1e6:.3f} km2')
        for lo, hi in ((1, 3), (4, 10), (11, 30), (1, 10), (1, 30)):
            print(f'   band l {lo}-{hi}: r {band_corr(V, H, LMAX, lo, hi):+.2f}')
    np.savez(Path(__file__).with_name('earth_degree_corr.npz'), r_e3=deg_corr(V, H3, LMAX), vd2020_alm=V,
             convention=np.array([best[1], best[2]]))


if __name__ == '__main__':
    main()
