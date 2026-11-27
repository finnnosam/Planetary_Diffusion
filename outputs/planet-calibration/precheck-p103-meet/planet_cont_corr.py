"""Pass 103 probe: how the planet's dynamic field meets its continents, by degree, beside Earth's (earth_cont_corr.txt).
Per arm and seed, snapshots at the tags given: the whole air-loaded field (dynamic_topography_air_m, swells included)
and its placed part (flow_topography_alm_m) expanded to l 1-30; the label-free continent (crust + bulk sediment >=
20 km) as an indicator function. Prints, means over seeds and tags: the degree correlation r_l at l 1-10 and by band
(l 1-3 / 4-10 / 11-30) of the whole field, the placed part and the rest with the continent; the same by band over the
oceanic nodes only (the reading Earth's oceans constrain); the field's area mean on the continent and on the ocean by
band. Earth (VD2020 / Holdt 2022 against CRUST1.0's c20): l 1-3 -0.62 / -0.31 (oceans -0.69 / -0.42), l 4-10 -0.06 /
-0.09 (oceans -0.06 / +0.01), l 11-30 -0.02 / +0.03. Read, not a rule. Usage (repo root): python3 planet_cont_corr.py
<arm> ... [--tags t+0000.0,t-0020.0,t-0040.0]"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_reviews.crust_review import sediment_bulk_m

LM = 30
BANDS = ((1, 3), (4, 10), (11, 30))


def deg_corr(A, B, lmax):
    l, m = mf.alm_l(lmax), mf.alm_m(lmax); w = np.where(m == 0, 1., 2.)
    num = np.bincount(l, w*(A*np.conj(B)).real, lmax+1)
    den = np.sqrt(np.bincount(l, w*np.abs(A)**2, lmax+1)*np.bincount(l, w*np.abs(B)**2, lmax+1))
    return np.where(den > 0, num/np.where(den > 0, den, 1), np.nan)


def band_corr(A, B, lmax, lo, hi):
    l, m = mf.alm_l(lmax), mf.alm_m(lmax); w = np.where(m == 0, 1., 2.); s = (l >= lo) & (l <= hi)
    return float((w[s]*(A[s]*np.conj(B[s])).real).sum()/np.sqrt((w[s]*np.abs(A[s])**2).sum()*(w[s]*np.abs(B[s])**2).sum()))


def corr(x, y, w):
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.average((x-xm)*(y-ym), weights=w)/np.sqrt(np.average((x-xm)**2, weights=w)*np.average((y-ym)**2, weights=w)))


def trim(A, lfrom, lto):
    """alm of degree <= lto from an array of lmax lfrom (the flow's 31 -> 30)."""
    out = np.zeros((lto+1)*(lto+2)//2, np.complex128)
    for m in range(lto+1):
        for l in range(m, lto+1):
            out[mf.alm_index(l, m, lto)] = A[mf.alm_index(l, m, lfrom)]
    return out


def main(argv):
    tags = ['t+0000.0']
    if '--tags' in argv:
        i = argv.index('--tags'); tags = argv[i+1].split(','); argv = argv[:i]+argv[i+2:]
    g = CubeSphere(217); a = g.area; L = mf.alm_l(LM)
    for arm in map(Path, argv):
        RL, RB, RO, MN = [], [], [], []
        seeds = sorted(int(p.name[4:]) for p in arm.glob('seed*') if (p/'snapshots').is_dir())
        for s in seeds:
            for t in tags:
                q = np.load(arm/f'seed{s}'/'snapshots'/f'{t}.npz')
                c = (q['crust_thickness_m']+sediment_bulk_m(q)) >= 20e3
                C = mf.point_alm(g.xyz, (c-c.mean()).astype(float)*a/g.radius**2, LM)
                W = mf.point_alm(g.xyz, q['dynamic_topography_air_m'].astype(float)*a/g.radius**2, LM)
                lm = int(q['flow_lmax']); P = trim(q['flow_topography_alm_m'].astype(np.complex128), lm, LM)
                Rr = W-P
                RL.append([deg_corr(F, C, LM)[1:11] for F in (W, P, Rr)])
                RB.append([[band_corr(F, C, LM, lo, hi) for lo, hi in BANDS] for F in (W, P, Rr)])
                o = ~c; ro, mn = [], []
                for lo, hi in BANDS:
                    s_ = (L >= lo) & (L <= hi)
                    f = mf.synthesis(np.where(s_, W, 0), LM, g.xyz); cc = mf.synthesis(np.where(s_, C, 0), LM, g.xyz)
                    ro.append(corr(f[o], cc[o], a[o]))
                    mn.append([(a[c]*f[c]).sum()/a[c].sum(), (a[o]*f[o]).sum()/a[o].sum()])
                RO.append(ro); MN.append(mn)
        RL, RB, RO, MN = map(np.array, (RL, RB, RO, MN))
        print(f'{arm.name} ({len(seeds)} seeds x {len(tags)} snapshots; means; se over runs in brackets)')
        for k, name in enumerate(('whole field', 'placed part', 'rest')):
            print(f'   {name:11} r_l l 1-10: ' + ' '.join(f'{v:+.2f}' for v in RL[:, k].mean(0)))
            print(f'   {name:11} bands l 1-3 / 4-10 / 11-30: ' + ' / '.join(
                f'{RB[:, k, j].mean():+.2f} ({RB[:, k, j].std()/np.sqrt(len(RB)):.2f})' for j in range(3)))
        print('   whole field, oceans only, bands: ' + ' / '.join(f'{RO[:, j].mean():+.2f}' for j in range(3)))
        print('   whole field area mean on continent / ocean (m) by band: ' + '; '.join(
            f'l {lo}-{hi} {MN[:, j, 0].mean():+.0f} / {MN[:, j, 1].mean():+.0f}' for j, (lo, hi) in enumerate(BANDS)))
        print('   by seed, whole field l 1-3 band: ' + ', '.join(f'{v:+.2f}' for v in RB[::len(tags), 0, 0]))


if __name__ == '__main__':
    main(sys.argv[1:])
