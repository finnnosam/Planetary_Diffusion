"""Pass 82 step-back, probe 4: static counterfactuals. On the planet's own 0 Myr states one thing is replaced by Earth's
(or removed), the columns are re-floated under water where the shore moves (a continental column under the sea sinks by
rho_w / (rho_m - rho_w) of its depth, as surface.elevation) and the sea is re-solved for the model's water volume.
No trajectory response (erosion, deposition, plates): this measures how much height each thing holds today, not what a
run with the change would give. Six seeds; mean (sd between seeds) of the change against the state as run.
  floor: every ocean-floor node raised by Earth's minus the planet's mean (sediment + crust terms + what the law leaves)
         in its age class (floor.law: Earth's residual only; floor.sed: the sediment term only; floor.crust)
  dyn0: the dynamic topography removed everywhere
  LT: continental L shifted so its area mean in each crust bin is Earth's seismic mean (Hoggard SL2013sv), spread kept
  LTfull: L a function of crust alone (Earth's bin means), the age spread removed
  youngthick: young (< 250 Myr) lithosphere under crust > 45 km set to Earth's 138 (45-55 km) / 197 km (> 55 km)
  cratonmean: every seed's continental mean L set to the frame's reference (130.68 km), a uniform shift"""
import sys, numpy as np, common as C
from planet_tectonics.isostasy import gdh1_depth, sea_level_for_volume
from planet_tectonics import lithosphere as lt
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM; tag = sys.argv[2] if len(sys.argv) > 2 else 't+0000.0'
AG = [(0, 10), (10, 30), (30, 60), (60, 100), (100, 150), (150, 400)]
TB = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]; TC = np.array([22.5, 27.5, 32.5, 37.5, 42.5, 50., 62.])
e = C.earth(); p = e['p']; rc, rm, rw = p['rho_crust'], p['rho_mantle'], p['rho_water']; V = p['water_volume_km3']*1e9; K = lt.LAMB_K_M_PER_KM
ds = np.load('outputs/planet-reference/earth-seafloor-age-n217.npz')['dated_share']
fl = (e['T'] < 20e3) & (ds >= .5) & (e['z'] < 0); oae = e['oa']
h = e['gsed']; phi0, c = p['phi0_ocean'], p['compaction_c_ocean_per_km']/1e3
sed_e = h-(p['rho_grain']-rw)*(h-phi0/c*(1-np.exp(-c*h)))/(rm-rw); crust_e = (e['Tc']-p['ocean_thick_m'])*(rm-rc)/(rm-rw)
res_e = e['z']-(-gdh1_depth(oae)+crust_e+sed_e)
EK = np.array([[C.wmean(t[m], e['a'][m]) for t in (sed_e, crust_e, res_e)] for m in [fl & (oae >= a) & (oae < b) for a, b in AG]])
ce = e['cont'] & ~e['ice']
LE = np.array([C.wmean(e['L'][ce & (e['tot'] >= a) & (e['tot'] < b)], e['a'][ce & (e['tot'] >= a) & (e['tot'] < b)]) for a, b in TB])
print('Earth floor terms by age class (sediment, crust, residual; m):', np.round(EK).tolist()); print('Earth seismic L by crust bin (km):', np.round(LE).tolist())
def resolve(P, zd_new):
    """Hydro-isostatic re-solve: dry-equivalent continental columns, re-floated under the new sea."""
    a, cont, s0 = P['a'], P['cont'], P['sea']; f = rm/(rm-rw)
    dry = np.where(cont & (P['zd'] < s0), s0+(P['zd']-s0)/f, P['zd'])+(zd_new-P['zd'])
    s = s0
    for _ in range(30):
        z = np.where(cont & (dry < s), s+(dry-s)*f, dry); s1 = sea_level_for_volume(z, a, V)
        if abs(s1-s) < .01: break
        s = s1
    return z, s1
def read(P, z, s):
    a = P['a']; l = z > s; q = C.wq(z[l]-s, a[l], [.25, .5, .9])
    return [s, a[l].sum()/a.sum(), *q]
SC = ['as run', 'floor', 'floor.law', 'floor.sed', 'floor.crust', 'dyn0', 'LT', 'LTfull', 'youngthick', 'cratonmean', 'floor+dyn0', 'floor+dyn0+LT', 'floor+dyn0+LT+youngthick']
R = {k: [] for k in SC}
for sd in C.SEEDS:
    P = C.planet(sd, tag, arm); a = P['a']; cont = P['cont']; oc = ~cont; zd = P['zd']
    cls = [oc & (P['oa'] >= lo) & (P['oa'] < hi) for lo, hi in AG]
    PK = np.array([[C.wmean(t[m], a[m]) if m.any() else 0. for t in (P['sedo'], P['crust_o'])] for m in cls])
    d_sed = np.zeros(zd.size); d_cr = np.zeros(zd.size); d_law = np.zeros(zd.size)
    for k, m in enumerate(cls):
        d_sed[m] = EK[k, 0]-PK[k, 0]; d_cr[m] = EK[k, 1]-PK[k, 1]; d_law[m] = EK[k, 2]
    cc = cont & (P['tot'] >= 20); L = P['L']
    LP = np.array([C.wmean(L[cc & (P['tot'] >= lo) & (P['tot'] < hi)], a[cc & (P['tot'] >= lo) & (P['tot'] < hi)]) if (cc & (P['tot'] >= lo) & (P['tot'] < hi)).any() else np.nan for lo, hi in TB])
    ok = np.isfinite(LP); shift = np.interp(P['tot'], TC[ok], (LE-LP)[ok])
    d_LT = np.where(cc, -K*shift, 0.); d_LTfull = np.where(cc, -K*(np.interp(P['tot'], TC, LE)-L), 0.)
    yt = cc & (P['age'] < 250) & (P['tot'] > 45); d_yt = np.where(yt, -K*(np.where(P['tot'] > 55, 197., 138.)-L), 0.); d_yt = np.minimum(d_yt, 0.)
    d_cm = np.where(cont, -K*(lt.LAMB_REF_KM-C.wmean(L[cont], a[cont])), 0.)
    mods = {'as run': 0., 'floor': d_sed+d_cr+d_law, 'floor.law': d_law, 'floor.sed': d_sed, 'floor.crust': d_cr, 'dyn0': -P['dyn'], 'LT': d_LT,
            'LTfull': d_LTfull, 'youngthick': d_yt, 'cratonmean': d_cm, 'floor+dyn0': d_sed+d_cr+d_law-P['dyn'],
            'floor+dyn0+LT': d_sed+d_cr+d_law-P['dyn']+d_LT, 'floor+dyn0+LT+youngthick': d_sed+d_cr+d_law-P['dyn']+d_LT+d_yt}
    for k in SC:
        z, s = resolve(P, zd+mods[k]); R[k].append(read(P, z, s))
    if sd == 1001: print(f'check seed 1001: re-solved sea as run {R["as run"][0][0]:.1f} m against the snapshot {P["sea"]:.1f}; planet L by crust bin {np.round(LP).tolist()}')
base = np.array(R['as run'])
print(f'\n{arm.split("/")[-1]} {tag}; Earth: land 0.277, p25 / p50 / p90 165 / 380 / 1414 m')
print('  scenario                  | sea on datum | land share | land p25 | p50 | p90 || change of p50 (sd seeds) | of p25 | of p90 | of the sea')
for k in SC:
    r = np.array(R[k]); m = r.mean(0); d = r-base
    print(f'  {k:25} | {m[0]:6.0f} | {m[1]:.3f} | {m[2]:5.0f} | {m[3]:5.0f} | {m[4]:5.0f} || {d[:, 3].mean():+5.0f} ({d[:, 3].std(ddof=1):3.0f}) | {d[:, 2].mean():+5.0f} | {d[:, 4].mean():+5.0f} | {d[:, 0].mean():+5.0f}')
