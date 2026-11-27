"""Pass 82 step-back, probe 6: journal 98's column_split re-read on today's arm, to see what its three shares depend
on. Same construction (sea to the datum; planet nodes reweighted to Earth's crust-thickness distribution; the rest is
'the column'), area-weighted, with (a) both orders of the two steps, (b) the thickness bins made finer (its bins;
2.5 km; 1 km), (c) total (crust + sediment) or crystalline thickness. Earth's distribution: CRUST1.0 on the nodes,
non-oceanic types, ice out (the domain of earth-crust-elevation-n217.json). Six seeds, mean of per-seed readings."""
import sys, numpy as np, common as C
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM
e = C.earth(); ae = e['a']; dom = ~e['oceanic'] & ~e['ice']
le = e['z'] > 0; E50 = C.wq(e['z'][le], ae[le], .5)
BINS = {'journal 98 bins': [0, 15, 20, 25, 30, 35, 40, 45, 55, 80], '2.5 km': [0, 15]+list(np.arange(17.5, 60.1, 2.5))+[80], '1 km': [0, 15]+list(np.arange(16, 60.1, 1.))+[80]}
print(f'Earth land p50 {E50:.0f} m')
for tname, key in (('total (crust + sediment)', 'tot'), ('crystalline', 'cry')):
    for bname, B in BINS.items():
        B = np.array(B, float); out = []
        TE = e['T']/1e3 if key == 'tot' else e['Tc']/1e3
        fe = np.array([ae[dom & (TE >= a) & (TE < b)].sum() for a, b in zip(B[:-1], B[1:])])/ae.sum()
        for s in C.SEEDS:
            P = C.planet(s, 't+0000.0', arm); a = P['a']; z = P['zd']; sea = P['sea']
            T = (P['T']+P['sed'])/1e3 if key == 'tot' else P['T']/1e3
            w = np.zeros(z.size)
            for (lo, hi), f in zip(zip(B[:-1], B[1:]), fe):
                m = (T >= lo) & (T < hi)
                if m.any(): w[m] = f/(a[m].sum()/a.sum())
            r = []
            for s0, wt in ((sea, None), (0., None), (sea, w), (0., w)):
                l = z > s0; ww = a[l]*(1 if wt is None else wt[l]); r.append(C.wq(z[l]-s0, ww, .5))
            out.append(r)
        r = np.mean(out, 0); sd = np.std(out, 0, ddof=1)
        print(f'{tname:26} | {bname:15} | as run {r[0]:.0f}, sea at the datum {r[1]:.0f}, Earth weights {r[2]:.0f}, both {r[3]:.0f} (sd {sd[3]:.0f}) | '
              f'sea first: sea {r[0]-r[1]:.0f}, crust {r[1]-r[3]:.0f}; crust first: crust {r[0]-r[2]:.0f}, sea {r[2]-r[3]:.0f}; column {r[3]-E50:.0f} | '
              f'shares of {r[0]-E50:.0f}: sea {100*(r[0]-r[1]+r[2]-r[3])/2/(r[0]-E50):.0f}%, crust {100*(r[1]-r[3]+r[0]-r[2])/2/(r[0]-E50):.0f}%, column {100*(r[3]-E50)/(r[0]-E50):.0f}%')
