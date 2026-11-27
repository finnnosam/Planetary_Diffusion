"""Pass 82 step-back, probe 7: the surface by crust class (label-free total thickness, crust + sediment): share of the
surface, mean elevation on the datum, and the room each holds below the datum; planet (six seeds, 0 and -250 Myr)
beside Earth (CRUST1.0 total on the nodes, ETOPO; oceanic-type share of each class printed). And the variance of land
elevation by term: cov(term, z) / var(z) over land (what high land has that low land has not)."""
import sys, numpy as np, common as C
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM
B = [0, 10, 15, 20, 25, 30, 35, 40, 45, 55, 100]
e = C.earth(); ae = e['a']/e['a'].sum(); ze = e['z']; TE = e['T']/1e3
print('class km | Earth: share / mean z / room / oceanic-type share | planet 0 Myr: share (sd) / mean z on datum / room | planet -250 Myr: share / mean z / room | room excess 0 Myr')
rows = {}
for tag in ('t+0000.0', 't-0250.0'):
    Ps = [C.planet(s, tag, arm) for s in C.SEEDS]; rows[tag] = []
    for lo, hi in zip(B[:-1], B[1:]):
        r = []
        for P in Ps:
            a = P['a']/P['a'].sum(); m = (P['tot'] >= lo) & (P['tot'] < hi)
            r.append([a[m].sum(), C.wmean(P['zd'][m], a[m]) if m.any() else np.nan, np.sum(a[m]*np.maximum(-P['zd'][m], 0))])
        rows[tag].append((np.nanmean(r, 0), np.nanstd(r, 0, ddof=1)))
tot = 0
for k, (lo, hi) in enumerate(zip(B[:-1], B[1:])):
    m = (TE >= lo) & (TE < hi); re = [ae[m].sum(), C.wmean(ze[m], ae[m]), np.sum(ae[m]*np.maximum(-ze[m], 0)), ae[m & e['oceanic']].sum()/ae[m].sum()]
    p0, s0 = rows['t+0000.0'][k]; p1, _ = rows['t-0250.0'][k]; tot += p0[2]-re[2]
    print(f'  {lo:3d}-{hi:3d} | {re[0]:.3f} / {re[1]:6.0f} / {re[2]:5.0f} / {re[3]:.2f} | {p0[0]:.3f} ({s0[0]:.3f}) / {p0[1]:6.0f} / {p0[2]:5.0f} | {p1[0]:.3f} / {p1[1]:6.0f} / {p1[2]:5.0f} | {p0[2]-re[2]:5.0f}')
print(f'  room excess summed {tot:.0f} m')
# variance of land elevation by term
cell = np.clip((90-e['lat']).astype(int), 0, 179)*360+np.clip((e['lon']+180).astype(int), 0, 359)
lf = np.bincount(cell, weights=(ze > 0)*e['a'], minlength=64800)/np.maximum(np.bincount(cell, weights=e['a'], minlength=64800), 1e-9)
def cov_shares(z, a, m, terms):
    a = a[m]/a[m].sum(); zc = z[m]-np.sum(a*z[m]); v = np.sum(a*zc*zc)
    return np.sqrt(v), [np.sum(a*zc*(t[m]-np.sum(a*t[m])))/v for t in terms]
land = (ze > 0) & ~e['ice'] & (lf[cell] >= .9)
print('\nland elevation: sd (m) and the share of its variance carried by each term, cov(term, z) / var(z): crust (Airy) | densification | sediment | lithosphere C | dyn | rest (Earth: residual)')
for cn in ('c_seis', 'c_age'):
    res = ze-(e['airy']+e['dens']+e['sedc']+e[cn]); sd, sh = cov_shares(ze, e['a'], land & np.isfinite(e[cn]), [e['airy'], e['dens'], e['sedc'], e[cn], 0*ze, res])
    print(f'  Earth (full-land cells, C from {cn}): sd {sd:.0f} | ' + ' | '.join(f'{x:+.2f}' for x in sh))
for tag in ('t+0000.0', 't-0250.0'):
    for F, nm in ((1., 'land as run'), (0., 'land above the datum')):
        o = []
        for s in C.SEEDS:
            P = C.planet(s, tag, arm); zs = P['zd']-F*P['sea']; m = (zs > 0) & P['cont']
            sd, sh = cov_shares(zs, P['a'], m, [P['airy'], P['dens'], P['sedc'], P['cterm'], P['dyn'], P['rest']]); o.append([sd]+sh)
        o = np.array(o); mm, ss = o.mean(0), o.std(0, ddof=1)
        print(f'  planet {tag}, {nm}: sd {mm[0]:.0f} ({ss[0]:.0f}) | ' + ' | '.join(f'{x:+.2f} ({y:.2f})' for x, y in zip(mm[1:], ss[1:])))
