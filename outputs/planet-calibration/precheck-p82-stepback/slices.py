"""Pass 82 step-back, probe 1: what the land's height is made of, term by term, at the same place in the land's own
distribution. For land nodes in a quantile slice of land elevation (p20-30, p45-55, p85-95, and all land), the area
mean of every term of the model's elevation law; planet (six seeds, 0 and -250 Myr) beside Earth (the same law on
CRUST1.0 / Hoggard inputs; what it leaves is Earth's residual). Each row sums exactly to the slice's mean elevation
above the sea, so planet minus Earth by column is an exact split of the excess at that quantile. Conditioning is on
elevation (ETOPO), never on CRUST1.0 thickness, so CRUST1.0's errors do not bias the means (no regression dilution).
Usage: python slices.py [arm [F]]"""
import sys, numpy as np, common as C
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM
F = float(sys.argv[2]) if len(sys.argv) > 2 else 1.   # 1: the sea as run; 0: the sea at the datum (land = above the datum)
SL = [('p20-30', .20, .30), ('p45-55', .45, .55), ('p85-95', .85, .95), ('all land', 0., 1.)]
COLS = ['z', 'airy', 'dens', 'sed', 'C', 'dyn', 'rest', '-sea', 'Tc_km', 'sed_m', 'L_km']
def rows(z, a, land, terms):
    out = {}
    lo_hi = {n: C.wq(z[land], a[land], [lo, hi]) for n, lo, hi in SL}
    for n, lo, hi in SL:
        q0, q1 = lo_hi[n]; m = land & (z >= q0) & (z <= q1) if n != 'all land' else land
        out[n] = [C.wmean(t[m], a[m]) for t in terms]
    return out
def show(name, R, sd=None):
    print(f'\n{name}\n  slice    | ' + ' | '.join(f'{c:>7}' for c in COLS))
    for n, _, _ in SL:
        print(f'  {n:8} | ' + ' | '.join(f'{v:7.1f}' if c.endswith('km') else f'{v:7.0f}' for c, v in zip(COLS, R[n])))
        if sd: print(f'  {"sd seeds":8} | ' + ' | '.join(f'{v:7.1f}' if c.endswith('km') else f'{v:7.0f}' for c, v in zip(COLS, sd[n])))
e = C.earth(); a = e['a']; z = e['z']
cell = np.clip((90-e['lat']).astype(int), 0, 179)*360+np.clip((e['lon']+180).astype(int), 0, 359)
lf = np.bincount(cell, weights=(z > 0)*a, minlength=64800)/np.maximum(np.bincount(cell, weights=a, minlength=64800), 1e-9)
ER = {}
for cname in ('c_age', 'c_seis', 'c_seis_clip'):
    for dom, land in (('all land, ice out', (z > 0) & ~e['ice']), ('full-land 1-degree cells', (z > 0) & ~e['ice'] & (lf[cell] >= .9)),
                      ('full-land cells, continental type, >= 20 km', (z > 0) & ~e['ice'] & (lf[cell] >= .9) & e['cont'])):
        ok = land & np.isfinite(e[cname])
        res = z-(e['airy']+e['dens']+e['sedc']+e[cname])
        R = rows(z, a, ok, [z, e['airy'], e['dens'], e['sedc'], e[cname], 0*z, res, 0*z, e['Tc']/1e3, e['S'], e['L']])
        ER[(cname, dom)] = R
        show(f'Earth, C from {cname}; {dom} (land share of surface {C.wmean(ok.astype(float), a):.3f}); "rest" = Earth\'s residual', R)
E0 = ER[('c_seis', 'full-land 1-degree cells')]
for tag in ('t+0000.0', 't-0250.0'):
    per = []
    for s in C.SEEDS:
        P = C.planet(s, tag, arm); zs = P['zd']-F*P['sea']; land = (zs > 0) & P['cont']
        per.append(rows(zs, P['a'], land, [zs, P['airy'], P['dens'], P['sedc'], P['cterm'], P['dyn'], P['rest'], 0*zs-F*P['sea'],
                                           P['T']/1e3, P['sh'], P['L']]))
    M = {n: np.mean([r[n] for r in per], 0) for n, _, _ in SL}; S = {n: np.std([r[n] for r in per], 0, ddof=1) for n, _, _ in SL}
    show(f'Planet {arm.split("/")[-1]} {tag}, sea at {F:g} x its level on the datum (six seeds: mean, sd between seeds)', M, S)
    print('  planet minus Earth (C from seismic L, full-land cells):')
    for n, _, _ in SL:
        print(f'  {n:8} | ' + ' | '.join(f'{v:7.1f}' if c.endswith('km') else f'{v:7.0f}' for c, v in zip(COLS, np.array(M[n])-np.array(E0[n]))))
