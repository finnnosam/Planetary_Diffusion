"""Pass 82 step-back, probe 10: seed by seed (0 Myr): the sea, the land, the continent's mean height on the datum, the
lithosphere under the land; land p50 against the sea across the six seeds and across the saved snapshots of each seed
(-40 ... 0 Myr). Earth beside it (continent >= 30 km: CRUST1.0 total, ice out)."""
import glob, numpy as np, common as C
e = C.earth(); ae = e['a']; m = (e['T'] >= 30e3) & ~e['ice']; l = (e['z'] > 0) & ~e['ice']
print(f"Earth: crust >= 30 km (ice out) share {ae[m].sum()/ae.sum():.3f}, mean elevation {C.wmean(e['z'][m], ae[m]):.0f} m, median {C.wq(e['z'][m], ae[m], .5):.0f}; "
      f"mean crust {C.wmean(e['T'][m], ae[m])/1e3:.1f} km; L under land (seismic) {C.wmean(e['L'][l], ae[l]):.1f} km")
print('seed | sea | land share | land p50 | crust >= 30 km: share, mean z on the datum, median, mean total km | L under land km | crust > 45 km share of continent | dyn mean on land')
R = []
for s in C.SEEDS:
    P = C.planet(s); a = P['a']; zs = P['zd']-P['sea']; land = zs > 0; k = P['tot'] >= 30; c20 = P['tot'] >= 20
    r = [P['sea'], a[land].sum()/a.sum(), C.wq(zs[land], a[land], .5), a[k].sum()/a.sum(), C.wmean(P['zd'][k], a[k]), C.wq(P['zd'][k], a[k], .5), C.wmean(P['tot'][k], a[k]),
         C.wmean(P['L'][land & P['cont']], a[land & P['cont']]), a[P['tot'] > 45].sum()/a[c20].sum(), C.wmean(P['dyn'][land], a[land])]
    R.append(r); print(s, ' | '.join(f'{v:.3f}' if abs(v) < 1 else f'{v:.0f}' if abs(v) > 200 else f'{v:.1f}' for v in r))
R = np.array(R); print('mean', np.round(R.mean(0), 3).tolist()); print('sd  ', np.round(R.std(0, ddof=1), 3).tolist())
b = np.polyfit(R[:, 0], R[:, 2], 1); print(f'across seeds at 0 Myr: land p50 = {b[1]:.0f} + {b[0]:.2f} x sea; corr {np.corrcoef(R[:, 0], R[:, 2])[0, 1]:.2f}; corr(p50, L under land) {np.corrcoef(R[:, 7], R[:, 2])[0, 1]:.2f}; corr(p50, crust > 45) {np.corrcoef(R[:, 8], R[:, 2])[0, 1]:.2f}; corr(p50, mean z of crust >= 30 km on the datum) {np.corrcoef(R[:, 4], R[:, 2])[0, 1]:.2f}')
X = np.c_[np.ones(6), R[:, 0], R[:, 4]]; co = np.linalg.lstsq(X, R[:, 2], rcond=None)[0]; fit = X@co
print(f'two-term fit: land p50 = {co[0]:.0f} + {co[1]:.2f} x sea + {co[2]:.2f} x (mean z of crust >= 30 km on the datum); rms misfit {np.sqrt(np.mean((fit-R[:, 2])**2)):.0f} m of sd {R[:, 2].std():.0f}')
sl = []
for s in C.SEEDS:
    pts = []
    for f in sorted(glob.glob(f'{C.ARM}/seed{s}/snapshots/t*.npz')):
        tag = f.split('/')[-1][:-4]
        if float(tag[1:]) < -40: continue
        P = C.planet(s, tag); zs = P['zd']-P['sea']; land = zs > 0; pts.append([P['sea'], C.wq(zs[land], P['a'][land], .5)])
    pts = np.array(pts); sl.append(np.polyfit(pts[:, 0], pts[:, 1], 1)[0] if np.ptp(pts[:, 0]) > 5 else np.nan)
print('within a seed over -40 ... 0 Myr (9 snapshots): slope of land p50 on the sea, by seed', np.round(sl, 2).tolist())
