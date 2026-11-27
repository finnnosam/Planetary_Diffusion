"""Pass 82 step-back, probe 11: the continents on the datum, with no sea and no land mask in the reading. Domain: total
crust (crust + sediment) >= 30 km (Earth: CRUST1.0, ice out; 0.278 of the surface; planet 0.329). Quantiles of
elevation on the datum, and in quantile slices of it the area mean of every term of the law (as slices.py; columns
under water are in: their water load is in 'rest'). Planet six seeds (0, -250 Myr) beside Earth."""
import sys, numpy as np, common as C
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM
Q = [.05, .10, .25, .50, .75, .90, .95]
SL = [('p10-20', .10, .20), ('p25-35', .25, .35), ('p45-55', .45, .55), ('p65-75', .65, .75), ('p85-95', .85, .95)]
COLS = ['z', 'airy', 'dens', 'sed', 'C', 'dyn', 'rest', 'tot_km', 'L_km']
def rows(z, a, dom, terms):
    out = {}
    for n, lo, hi in SL:
        q0, q1 = C.wq(z[dom], a[dom], [lo, hi]); m = dom & (z >= q0) & (z <= q1); out[n] = [C.wmean(t[m], a[m]) for t in terms]
    return out
def show(name, R):
    print(f'{name}\n  slice    | ' + ' | '.join(f'{c:>7}' for c in COLS))
    for n, _, _ in SL: print(f'  {n:8} | ' + ' | '.join(f'{v:7.1f}' if c.endswith('km') else f'{v:7.0f}' for c, v in zip(COLS, R[n])))
e = C.earth(); a = e['a']; z = e['z']; dom = (e['T'] >= 30e3) & ~e['ice'] & np.isfinite(e['c_seis'])
qe = C.wq(z[dom], a[dom], Q); print('quantiles of elevation on the datum, crust >= 30 km: ' + ' / '.join(f'p{int(100*q)}' for q in Q))
print('  Earth  : ' + ' / '.join(f'{v:.0f}' for v in qe) + f'   (share under the datum {a[dom & (z < 0)].sum()/a[dom].sum():.3f})')
res = z-(e['airy']+e['dens']+e['sedc']+e['c_seis'])
E = rows(z, a, dom, [z, e['airy'], e['dens'], e['sedc'], e['c_seis'], 0*z, res, e['tot'], e['L']])
for tag in ('t+0000.0', 't-0250.0'):
    qs, per, sub = [], [], []
    for s in C.SEEDS:
        P = C.planet(s, tag, arm); d = (P['tot'] >= 30) & P['cont']; zd = P['zd']
        qs.append(C.wq(zd[d], P['a'][d], Q)); sub.append(P['a'][d & (zd < 0)].sum()/P['a'][d].sum())
        per.append(rows(zd, P['a'], d, [zd, P['airy'], P['dens'], P['sedc'], P['cterm'], P['dyn'], P['rest'], P['tot'], P['L']]))
    qs = np.array(qs); print(f'  planet {tag}: ' + ' / '.join(f'{v:.0f}' for v in qs.mean(0)) + f'   (share under the datum {np.mean(sub):.3f}); sd between seeds ' + ' / '.join(f'{v:.0f}' for v in qs.std(0, ddof=1)))
    print('  planet - Earth: ' + ' / '.join(f'{v:+.0f}' for v in qs.mean(0)-qe))
    M = {n: np.mean([r[n] for r in per], 0) for n, _, _ in SL}
    show(f'\nplanet {tag} (six-seed mean)', M)
    if tag == 't+0000.0': show('Earth (C from seismic L; rest = Earth\'s residual)', E)
    print('  planet minus Earth:')
    for n, _, _ in SL: print(f'  {n:8} | ' + ' | '.join(f'{v:7.1f}' if c.endswith('km') else f'{v:7.0f}' for c, v in zip(COLS, np.array(M[n])-np.array(E[n]))))
