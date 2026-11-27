"""Pass 82 step-back, probe 2: the two hypsometric curves on one datum. Earth: ETOPO (sea = datum). Planet: elevation +
sea level (the datum the model's columns are referenced on: 38 km of crust at +300 m, ridges at -2,600 m).
(a) elevation at the same rank from the top of the surface; (b) share of the surface per elevation band about each
planet's own sea, and about the datum; (c) land p25/p50/p90 and land share as the planet's sea is moved from where it
is to the datum (nothing else changed). Six seeds pooled per node weight; per-seed spread printed."""
import sys, numpy as np, common as C
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM
e = C.earth(); ae = e['a']/e['a'].sum(); ze = e['z']
RANK = [.01, .03, .05, .08, .10, .1385, .1535, .20, .248, .277, .307, .34, .37, .40, .45]
def top(z, a, r): return C.wq(-z, a, r)*-1
for tag in ('t+0000.0', 't-0250.0'):
    Ps = [C.planet(s, tag, arm) for s in C.SEEDS]
    zs = [P['zd'] for P in Ps]; seas = np.array([P['sea'] for P in Ps]); a = Ps[0]['a']/Ps[0]['a'].sum()
    print(f'\n=== {tag}: sea on the datum {seas.mean():.0f} m (sd {seas.std(ddof=1):.0f}; by seed {np.round(seas)})')
    print('(a) elevation on the datum at the same rank from the top (share of the surface above it)')
    print('   rank | Earth | planet mean (sd seeds) | planet - Earth')
    for r in RANK:
        v = np.array([top(z, a, r) for z in zs]); ev = top(ze, ae, r)
        print(f'  {r:6.4f} | {ev:6.0f} | {v.mean():6.0f} ({v.std(ddof=1):4.0f}) | {v.mean()-ev:6.0f}')
    B = [-6000, -3000, -1000, -500, -200, 0, 200, 500, 1000, 2000, 9000]
    print('(b) share of the surface (%) per elevation band: Earth about its sea | planet about its sea | planet about the datum')
    for lo, hi in zip(B[:-1], B[1:]):
        fe = ae[(ze >= lo) & (ze < hi)].sum()
        fs = np.mean([a[(z-s >= lo) & (z-s < hi)].sum() for z, s in zip(zs, seas)]); fd = np.mean([a[(z >= lo) & (z < hi)].sum() for z in zs])
        print(f'  {lo:6d} ... {hi:5d} m | {100*fe:5.1f} | {100*fs:5.1f} | {100*fd:5.1f}')
    print('(c) the sea moved toward the datum, nothing else changed (sea level on the datum | land share | land p25 / p50 / p90 above it)')
    le = ze > 0; print(f'  Earth            | {ae[le].sum():.3f} | ' + ' / '.join(f'{v:.0f}' for v in C.wq(ze[le], ae[le], [.25, .5, .9])))
    for f in (1., .75, .5, .25, 0.):
        out = []
        for z, s in zip(zs, seas):
            l = z > f*s; out.append([a[l].sum(), *C.wq(z[l]-f*s, a[l], [.25, .5, .9])])
        o = np.mean(out, 0); sd = np.std(out, 0, ddof=1)
        print(f'  sea at {f*seas.mean():5.0f} m | {o[0]:.3f} | ' + ' / '.join(f'{v:.0f}' for v in o[1:]) + f'   (sd of p50 between seeds {sd[2]:.0f})')
