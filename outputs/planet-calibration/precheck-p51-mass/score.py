"""Pass 51: score the candidates of P51_RULE.txt on E3 (sh08): T1 rms l 1-31, T2 geoid VR l 2-8, the geoid's
preferred mass scale over l 4-9, E3 / tomographic proxy by bin (pass 49's ratios times the new / old mass per bin)."""
import json, glob
import numpy as np
from planet_tectonics import mantle_flow as mf
from planet_reviews import mantle_flow_earth as mfe
from planet_grid.paths import DATA

A = 6371e3; RHO_M = 3300.
E = {k: np.load(f'e3_area_{k}.npz') for k in (12, 15)}
RAD = E[12]['radii']; L = int(E[12]['lmax']); Z = (A-RAD)/1e3
base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08'); K = mf.FlowKernels(A, RAD, base, eta, L)
li, mi = mf.alm_l(L), mf.alm_m(L); W = np.where(mi == 0, 1., 2.)


def read_geoid(p):
    tok = open(p).read().split(); lm = int(tok[0]); v = np.array(tok[6:], float).reshape(-1, 2)
    return mfe.real_to_alm(v, lm)


OBS = {k: read_geoid(DATA/'mantle'/f'egm2008-hc-geoid.{f}.ab') for k, f in (('chambat', 'chambat.31'), ('nakiboglu', 'nakiboglu'))}


def vr(mod, obs, Ls):
    s = np.isin(li, Ls); return float(1-(W[s]*abs(obs[s]-mod[s])**2).sum()/(W[s]*abs(obs[s])**2).sum())


def scale(mod, obs, Ls):
    s = np.isin(li, Ls); return float((W[s]*(obs[s]*np.conj(mod[s])).real).sum()/(W[s]*abs(mod[s])**2).sum())


FD = mf.alpha_ratio(Z)
m1 = {f: json.load(open(f)) for f in glob.glob('m1_*.json')}
c = m1['m1_1600_63.json']; assert np.allclose(c['z_km'], Z)
F1c = 0.5*(np.array(c['F_noconduction'])+np.array(c['F_diffused']))
F1lo = np.min([np.minimum(v['F_noconduction'], v['F_diffused']) for v in m1.values()], 0)
F1hi = np.max([np.maximum(v['F_noconduction'], v['F_diffused']) for v in m1.values()], 0)
m2 = {f: np.array(json.load(open(f))['dm_kg_m2']) for f in glob.glob('m2_*.json')}
M2c = m2['m2_0.5_7.1_60.json']; M2lo = np.min(list(m2.values()), 0)
um = (Z < 670)
PROXY_RATIO = {'660-1000': 1.7, '1000-1500': 1.6, '1500-2000': 2.0, '2000-2500': 2.3, '2500-2850': 1.0}   # journal 70, with D's alpha


def run(name, F, M2=None, um_scale=1., sink=12):
    e = E[sink]; F = F*np.where(um, um_scale, 1.)
    loads = e['loads_0']*F[:, None]
    if M2 is not None:
        loads = loads+e['area_loads_0']*M2[:, None]
    H, _ = K.surface(loads); N = K.geoid(loads)
    P = mf.spectrum(H/RHO_M, L); rms = float(np.sqrt(P[1:].sum()))
    mb = e['mass_by_layer_0']*F+(0 if M2 is None else e['area_by_layer_0']*M2)
    mb0 = E[12]['mass_by_layer_0']*FD
    prox = {}
    for b, r in PROXY_RATIO.items():
        lo, hi = map(float, b.split('-')); s = (Z > lo) & (Z < hi); prox[b] = r*mb[s].sum()/mb0[s].sum()
    row = dict(rms_m=rms, share_l1_3=float(P[1:4].sum()/P[1:].sum()), T1=bool(150 <= rms <= 300),
               VR_l2_8={k: vr(N, o, range(2, 9)) for k, o in OBS.items()},
               s_geoid_l4_9=scale(N, OBS['chambat'], range(4, 10)), s_topo=220./rms, held_kg=float(mb.sum()), E3_over_proxy=prox)
    row['T2'] = bool(min(row['VR_l2_8'].values()) >= 0.5)
    print(f"{name:42s} rms {rms:5.0f} m (l1-3 {row['share_l1_3']:.2f}) T1 {'pass' if row['T1'] else 'FAIL'}; VR l2-8 "
          f"{row['VR_l2_8']['chambat']:+.2f}/{row['VR_l2_8']['nakiboglu']:+.2f} T2 {'pass' if row['T2'] else 'FAIL'}; geoid s l4-9 "
          f"{row['s_geoid_l4_9']:.2f}, topo s {row['s_topo']:.2f}; E3/proxy " + ' '.join(f'{v:.1f}' for v in prox.values()))
    return row


out = {}
out['D (alpha_ambient / alpha_GDH1)'] = run('D (alpha_ambient / alpha_GDH1)', FD)
out['M1 central'] = run('M1 central', F1c)
out['M1 low / high'] = [run('M1 low end', F1lo), run('M1 high end', F1hi)]
out['M2 central (on D)'] = run('M2 central (on D)', FD, M2c)
out['M1 + M2 central'] = run('M1 + M2 central', F1c, M2c)
out['M3 sink 15 mm/yr (on D)'] = run('M3 sink 15 mm/yr (on D)', FD, sink=15)
out['M1 + M2 central, sink 15'] = run('M1 + M2 central, sink 15', F1c, M2c, sink=15)
out['most favourable'] = run('most favourable (M1 low, M2 low, UM -10%, 15)', F1lo, M2lo, 0.9, 15)
json.dump(out, open('score.json', 'w'), indent=1)
print('factor by depth: z, D, M1 central, M1 low-high')
for j in range(0, len(Z), 5):
    print(f'  {Z[j]:5.0f} {FD[j]:.3f} {F1c[j]:.3f} {F1lo[j]:.3f}-{F1hi[j]:.3f}')
