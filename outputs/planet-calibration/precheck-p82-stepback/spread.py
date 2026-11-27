"""Pass 82 step-back, probe 14: what holds the spread of the continents' elevation on the datum (crust + sediment
>= 30 km; Earth p10 / p25 / p50 / p75 / p90 and p75 - p25 beside). The planet's 0 Myr columns with one term changed,
no re-floating (dry terms): dyn removed; L's crust-bin means set to Earth's (LT) or L a function of crust (LTfull); the
flexural part (rest) removed; crust above 45 km cut to 45 km in the Airy and densification terms (the thick tail's
height); sediment redistributed within the class so that its term follows Earth's by elevation rank (same class mean)."""
import numpy as np, common as C
from planet_tectonics import lithosphere as lt
from planet_tectonics.deep_crust import densified_mass
from planet_tectonics.isostasy import sediment_column
Q = [.10, .25, .50, .75, .90]; K = lt.LAMB_K_M_PER_KM
TB = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]; TC = np.array([22.5, 27.5, 32.5, 37.5, 42.5, 50., 62.])
e = C.earth(); ae = e['a']; ze = e['z']; de = (e['T'] >= 30e3) & ~e['ice']; qe = C.wq(ze[de], ae[de], Q)
ce = e['cont'] & ~e['ice']; LE = np.array([C.wmean(e['L'][ce & (e['tot'] >= a) & (e['tot'] < b)], ae[ce & (e['tot'] >= a) & (e['tot'] < b)]) for a, b in TB])
# Earth's sediment term by elevation rank within the class (for the redistribution)
RR = np.linspace(0, 1, 21); o = np.argsort(ze[de]); cw = np.cumsum(ae[de][o]); re = np.empty(de.sum()); re[o] = (cw-.5*ae[de][o])/cw[-1]
sed_by_rank = np.array([C.wmean(e['sedc'][de][(re >= a) & (re <= b)], ae[de][(re >= a) & (re <= b)]) for a, b in zip(RR[:-1], RR[1:])]); rc_ = (RR[:-1]+RR[1:])/2
print('quantiles on the datum p10 / p25 / p50 / p75 / p90 | p75 - p25 | p50 - p10')
print('  Earth                      | ' + ' / '.join(f'{v:6.0f}' for v in qe) + f' | {qe[3]-qe[1]:5.0f} | {qe[2]-qe[0]:5.0f}')
names = ['as run', 'dyn removed', 'LT', 'LTfull', 'flexure removed', 'crust cut at 45 km', 'sediment by rank as Earth', 'LT + dyn removed', 'LT + dyn + sediment by rank', 'all five']
out = {k: [] for k in names}
for s in C.SEEDS:
    P = C.planet(s); a = P['a']; z = P['zd']; p = P['p']; rm, rcr = p['rho_mantle'], p['rho_crust']; d = (P['tot'] >= 30) & P['cont']; L = P['L']
    LP = np.array([C.wmean(L[d & (P['tot'] >= lo) & (P['tot'] < hi)], a[d & (P['tot'] >= lo) & (P['tot'] < hi)]) if (d & (P['tot'] >= lo) & (P['tot'] < hi)).any() else np.nan for lo, hi in TB])
    cc = P['cont'] & (P['tot'] >= 20); LPa = np.array([C.wmean(L[cc & (P['tot'] >= lo) & (P['tot'] < hi)], a[cc & (P['tot'] >= lo) & (P['tot'] < hi)]) for lo, hi in TB])
    d_LT = -K*np.interp(P['tot'], TC, LE-LPa); d_LTf = -K*(np.interp(P['tot'], TC, LE)-np.nan_to_num(L))
    Tcut = np.minimum(P['T'], 45e3); sm_dry = sediment_column(P['sed'], P['cont'], p)[1]
    d_cut = (Tcut-P['T'])*(rm-rcr)/rm-(densified_mass(Tcut, sm_dry, P['cont'], p)-densified_mass(P['T'], sm_dry, P['cont'], p))/rm
    o = np.argsort(z[d]); cw = np.cumsum(a[d][o]); r = np.empty(d.sum()); r[o] = (cw-.5*a[d][o])/cw[-1]
    tgt = np.interp(r, rc_, sed_by_rank); tgt = tgt*C.wmean(P['sedc'][d], a[d])/C.wmean(tgt, a[d]); d_sed = np.zeros(z.size); d_sed[d] = tgt-P['sedc'][d]
    mods = {'as run': 0., 'dyn removed': -P['dyn'], 'LT': d_LT, 'LTfull': d_LTf, 'flexure removed': -P['rest'], 'crust cut at 45 km': d_cut, 'sediment by rank as Earth': d_sed,
            'LT + dyn removed': d_LT-P['dyn'], 'LT + dyn + sediment by rank': d_LT-P['dyn']+d_sed, 'all five': d_LT-P['dyn']+d_sed+d_cut-P['rest']}
    for k in names: out[k].append(C.wq((z+mods[k])[d], a[d], Q))
for k in names:
    q = np.mean(out[k], 0); sd = np.std([x[3]-x[1] for x in out[k]], ddof=1)
    print(f'  {k:26} | ' + ' / '.join(f'{v:6.0f}' for v in q) + f' | {q[3]-q[1]:5.0f} (sd {sd:3.0f}) | {q[2]-q[0]:5.0f}')
