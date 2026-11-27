"""Pass 50: D's viscosity profiles judged on Earth (P50_RULE.txt). E3's present-day record (e3_loads.npz) and E1's
tomographic loads through each candidate's kernels: residual topography (T1) and the geoid against EGM2008's
non-hydrostatic geoid, degrees 2-8 (T2). Writes profiles.json and prints the table."""
import json, sys
import numpy as np
from planet_tectonics import mantle_flow as mf
from planet_reviews import mantle_flow_earth as mfe
from planet_grid.paths import DATA

A = 6371e3; C = mf.CMB_RADIUS_M; RHO_M = 3300.
SEATREE = sys.argv[1] if len(sys.argv) > 1 else '/home/claude/seatree/python3/data/hc'


def layered(rows, um=1e21):
    """rows: (top depth km, relative viscosity) from the CMB up; returns (base fractions, eta) in the HC layout."""
    base = [C/A]+[(A-z*1e3)/A for z, _ in rows[1:]]
    return np.array(base), np.array([v*um for _, v in rows])


PROFILES = {
    'sh08': mf.read_viscosity(DATA/'mantle'/'visc.sh08'),
    'R93': layered([(2891, 37.), (670, 1.), (100, 11.)]),
    'LR98': layered([(2891, 50.), (670, 1.), (130, 10.)]),
    'H84': layered([(2891, 30.), (670, 1.)]),
}


def read_geoid(path):
    tok = open(path).read().split()
    lmax = int(tok[0]); v = np.array(tok[6:], float).reshape(-1, 2)
    assert v.shape[0] == (lmax+1)*(lmax+2)//2
    return lmax, mfe.real_to_alm(v, lmax)


LMAX = 31
li, mi = mf.alm_l(LMAX), mf.alm_m(LMAX)
W = np.where(mi == 0, 1., 2.)


def vr(mod, obs, L):
    s = np.isin(li, L)
    return float(1-(W[s]*np.abs(obs[s]-mod[s])**2).sum()/(W[s]*np.abs(obs[s])**2).sum())


def corr(x, y, L):
    s = np.isin(li, L)
    return float((W[s]*(x[s]*np.conj(y[s])).real).sum()/np.sqrt((W[s]*np.abs(x[s])**2).sum()*(W[s]*np.abs(y[s])**2).sum()))


def rms(x, L):
    s = np.isin(li, L)
    return float(np.sqrt((W[s]*np.abs(x[s])**2).sum()/(4*np.pi)))


OBS = {k: read_geoid(f'{SEATREE}/egm2008-hc-geoid.{f}.ab')[1] for k, f in (('chambat', 'chambat.31'), ('nakiboglu', 'nakiboglu'))}
E3 = np.load('e3_loads.npz'); RAD = E3['radii']; FAC = mf.alpha_ratio((A-RAD)/1e3)
L28 = range(2, 9)


def geoid_block(N):
    out = {'rms_l2_8_m': rms(N, L28)}
    for k, o in OBS.items():
        out[k] = {'VR_l2_8': vr(N, o, L28), 'corr_l2_8': corr(N, o, L28), 'amp_ratio_l2_8': rms(N, L28)/rms(o, L28),
                  'VR_l2_15': vr(N, o, range(2, 16)), 'corr_by_degree': {l: corr(N, o, [l]) for l in range(2, 13)},
                  'ratio_by_degree': {l: rms(N, [l])/rms(o, [l]) for l in range(2, 13)}}
    return out


def evaluate(base, eta, full=True):
    K = mf.FlowKernels(A, RAD, base, eta, LMAX)
    res = {}
    for t in ((0, 10, 20) if full else (0,)):
        loads = E3[f'loads_{t}']*FAC[:, None]
        H, V = K.surface(loads)
        P = mf.spectrum(H/RHO_M, LMAX)
        row = {'topo_rms_m': float(np.sqrt(P[1:].sum())), 'share_l1_3': float(P[1:4].sum()/P[1:].sum()),
               'flow_rms_cm_yr': float(np.sqrt(mf.gradient_power(V, LMAX))*3.15576e7*100)}
        row['geoid'] = geoid_block(K.geoid(loads))
        if t == 0 and full:
            for rt in (3300., 2270.):
                KN, KHs = K.geoid_kernels(True, rt)
                Ns = (KN[li, :].T*loads).sum(0); Hs = (KHs[li, :].T*loads).sum(0)
                Ps = mf.spectrum(Hs/RHO_M, LMAX)
                row[f'self_grav_{rt:.0f}'] = {'topo_rms_m': float(np.sqrt(Ps[1:].sum())), 'geoid': geoid_block(Ns)}
        res[t] = row
    if full:
        for key in (('smean2', True, 300.), ('smean2', False, 300.)):
            _, L = mfe.density_loads(*key)
            Kp = mf.FlowKernels(A, [A-z*1e3 for z, _ in L], base, eta, LMAX)
            ld = np.stack([a for _, a in L])
            H, _ = Kp.surface(ld); P = mf.spectrum(H/RHO_M, LMAX)
            res['proxy_fast' if key[1] else 'smean2_all'] = {'topo_rms_m': float(np.sqrt(P[1:].sum())),
                                                             'geoid': geoid_block(Kp.geoid(ld))}
    return res


def main():
    out = {'observed': {k: {'rms_l2_8_m': rms(o, L28), 'rms_l2_31_m': rms(o, range(2, 32)),
                            'deg_rms_m': {l: rms(o, [l]) for l in range(2, 9)}} for k, o in OBS.items()}}
    print('observed non-hydrostatic geoid:', json.dumps(out['observed']))
    for name, (base, eta) in PROFILES.items():
        r = evaluate(base, eta)
        out[name] = r
        g0 = r[0]['geoid']; T1 = .15e3 <= r[0]['topo_rms_m'] <= .30e3
        T2 = all(g0[k]['VR_l2_8'] >= .5 for k in OBS)
        print(f"\n{name}: E3 topo {r[0]['topo_rms_m']:.0f} m (l1-3 {r[0]['share_l1_3']:.2f}) T1 {'pass' if T1 else 'FAIL'}; geoid rms l2-8 "
              f"{g0['rms_l2_8_m']:.1f} m; VR l2-8 chambat {g0['chambat']['VR_l2_8']:+.2f} nakiboglu {g0['nakiboglu']['VR_l2_8']:+.2f} "
              f"T2 {'pass' if T2 else 'FAIL'}; corr {g0['chambat']['corr_l2_8']:+.2f}/{g0['nakiboglu']['corr_l2_8']:+.2f}, amp ratio "
              f"{g0['chambat']['amp_ratio_l2_8']:.2f}/{g0['nakiboglu']['amp_ratio_l2_8']:.2f}; VR l2-15 {g0['chambat']['VR_l2_15']:+.2f}; flow {r[0]['flow_rms_cm_yr']:.2f} cm/yr")
        print('   corr by degree (chambat):', ' '.join(f"{l}:{v:+.2f}" for l, v in g0['chambat']['corr_by_degree'].items()))
        print('   amp ratio by degree     :', ' '.join(f"{l}:{v:.2f}" for l, v in g0['chambat']['ratio_by_degree'].items()))
        for rt in ('3300', '2270'):
            s = r[0][f'self_grav_{rt}']
            print(f"   self-grav (surface {rt}): topo {s['topo_rms_m']:.0f} m, geoid VR l2-8 {s['geoid']['chambat']['VR_l2_8']:+.2f}/{s['geoid']['nakiboglu']['VR_l2_8']:+.2f}, rms {s['geoid']['rms_l2_8_m']:.1f} m")
        for t in (10, 20):
            print(f"   E3 at {t} Ma: topo {r[t]['topo_rms_m']:.0f} m, geoid VR l2-8 {r[t]['geoid']['chambat']['VR_l2_8']:+.2f}/{r[t]['geoid']['nakiboglu']['VR_l2_8']:+.2f}")
        for k in ('proxy_fast', 'smean2_all'):
            s = r[k]
            print(f"   {k}: topo {s['topo_rms_m']:.0f} m, geoid rms {s['geoid']['rms_l2_8_m']:.1f} m, VR l2-8 {s['geoid']['chambat']['VR_l2_8']:+.2f}/{s['geoid']['nakiboglu']['VR_l2_8']:+.2f}, corr {s['geoid']['chambat']['corr_l2_8']:+.2f}")
    json.dump(out, open('profiles.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
