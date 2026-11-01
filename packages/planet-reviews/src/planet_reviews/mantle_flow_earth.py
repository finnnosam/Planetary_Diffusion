"""D's Earth reference (journal 66, sections 3 and 5: E1, E2): the free-slip flow that Earth's own tomography drives
through planet_tectonics.mantle_flow (the planet's kernels, the same viscosity file), its dynamic-topography spectrum
and surface flow, and traction toward it in the torque law on Bird's plates. Reference numbers, not bands.

  python -m planet_reviews.mantle_flow_earth            E1 for the five density variants (a few seconds each)
  python -m planet_reviews.mantle_flow_earth --traction E1 and E2 (the traction table; Nelder-Mead refits, minutes)

Earth densities: S-wave tomography (SMEAN2, S40RTS; data/mantle/, from SEATREE) x dln rho / dln Vs = 0.2 x PREM-like
density. The slab proxy is the fast (positive density) part below the cut: below ~250 km keels, Earth's fast anomalies
are slabs. All of it below the cut also carries the hot upwellings."""
import argparse
import numpy as np

from planet_grid.paths import DATA
from planet_tectonics import mantle_flow as mf

YR_S = 3.15576e7
DLNRHO_DLNVS = .2
# a PREM-like mantle density with depth (km, kg/m^3), for converting dln rho to rho
PREM_Z = np.array([0, 220, 400, 670, 771, 1000, 1500, 2000, 2500, 2891.])
PREM_RHO = np.array([3380, 3360, 3540, 3990, 4440, 4580, 4800, 5050, 5300, 5560.])
MODELS = {'smean2': 'smean2.31.m.ab', 's40rts': 's40rts.31.m.ab'}
VARIANTS = (('smean2', True, 300.), ('s40rts', True, 300.), ('smean2', False, 300.), ('smean2', True, 100.),
            ('smean2', False, 100.))
# the pass-45 cloud pre-check (journal 66, section 3): rms (m, air-loaded), share l 1-3, P30/P2, surface flow (cm/yr)
PASS45 = {('smean2', True, 300.): (222., .84, 6e-4, 1.35), ('s40rts', True, 300.): (229., .82, 6e-4, 1.34),
          ('smean2', False, 300.): (417., .88, 3e-4, 2.58), ('smean2', True, 100.): (509., .48, 4e-3, 1.76),
          ('smean2', False, 100.): (850., .60, 2e-3, 3.34)}


def read_ab(path):
    """Becker & Boschi short format: number of layers; per layer depth (km), lmax, then (A, B) for l = 0..lmax,
    m = 0..l (fully normalised real harmonics). Returns lmax and [(depth km, (n, 2) array)] sorted by depth."""
    tok = open(path).read().split()
    nl = int(tok[0]); i = 1; out = []
    for _ in range(nl):
        z = float(tok[i]); lmax = int(tok[i+1]); i += 2; n = (lmax+1)*(lmax+2)//2
        v = np.array(tok[i:i+2*n], float).reshape(n, 2); i += 2*n
        out.append((z, v))
    out.sort(key=lambda t: t[0])
    return lmax, out


def real_to_alm(v, lmax):
    """Real (A, B) coefficients, ordered l-major, to ducc0's complex m-major alm (orthonormal)."""
    a = np.zeros((lmax+1)*(lmax+2)//2, complex); k = 0
    for l in range(lmax+1):
        for m in range(l+1):
            A, B = v[k]; k += 1
            a[mf.alm_index(l, m, lmax)] = A if m == 0 else (A-1j*B)/np.sqrt(2)
    return a


def density_loads(model, fast_only, z_cut_km):
    """[(depth km, alm of the layer's surface density, kg/m^2)] below the cut: dlnVs (%) x 0.2 x rho(z) x layer
    thickness (layers bounded halfway between the model's depths); fast_only keeps the positive part (on a
    Gauss-Legendre grid)."""
    import ducc0
    lmax, L = read_ab(DATA/'mantle'/MODELS[model])
    zs = np.array([z for z, _ in L]); edges = np.concatenate([[0], (zs[1:]+zs[:-1])/2, [2891.]])
    out = []
    for k, (z, v) in enumerate(L):
        if z < z_cut_km:
            continue
        a = real_to_alm(v, lmax)*DLNRHO_DLNVS*.01*np.interp(z, PREM_Z, PREM_RHO)
        if fast_only:
            m = ducc0.sht.synthesis_2d(alm=a[None], spin=0, lmax=lmax, geometry='GL', ntheta=2*lmax+2, nphi=4*lmax+4)[0]
            a = ducc0.sht.analysis_2d(map=np.maximum(m, 0)[None], spin=0, lmax=lmax, geometry='GL')[0]
        out.append((z, a*(edges[k+1]-max(edges[k], z_cut_km))*1e3))
    return lmax, out


def earth_flow(model='smean2', fast_only=True, z_cut_km=300., rho_m=3300.):
    """(lmax, H (m, air-loaded topography coefficients), V (m/s, velocity potential)) of the flow Earth's density
    drives, with the planet's kernels (radius 6371 km, data/mantle/visc.sh08)."""
    lmax, loads = density_loads(model, fast_only, z_cut_km)
    base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08')
    K = mf.FlowKernels(6371e3, [6371e3-z*1e3 for z, _ in loads], base, eta, lmax)
    H, V = K.surface(np.stack([a for _, a in loads]))
    return lmax, H/rho_m, V


def summary(lmax, H, V):
    P = mf.spectrum(H, lmax)
    return {'rms_m': float(np.sqrt(P[1:].sum())), 'share_l1_3': float(P[1:4].sum()/P[1:].sum()),
            'p30_p2': float(P[min(30, lmax)]/P[2]), 'deg_rms_m': {l: float(np.sqrt(P[l])) for l in (1, 2, 3, 4, 6, 8, 12, 20, 30)},
            'flow_rms_cm_yr': float(np.sqrt(mf.gradient_power(V, lmax))*YR_S*100)}


def e1():
    rows = []
    for key in VARIANTS:
        s = summary(*earth_flow(*key))
        rows.append((key, s))
    return rows


def e1_text(rows):
    out = ['E1: flow from Earth density through the planet kernels (air-loaded dynamic topography; pass-45 numbers in brackets)']
    for (model, fast, zc), s in rows:
        ref = PASS45[(model, fast, zc)]
        out.append(f"  {model} {'fast part' if fast else 'all'} > {zc:.0f} km: rms {s['rms_m']:.0f} m ({ref[0]:.0f}); share l1-3 "
                   f"{s['share_l1_3']:.2f} ({ref[1]:.2f}); P30/P2 {s['p30_p2']:.1e} ({ref[2]:.0e}); flow rms "
                   f"{s['flow_rms_cm_yr']:.2f} cm/yr ({ref[3]:.2f}); degree rms "
                   + ' '.join(f'{l}:{v:.0f}' for l, v in s['deg_rms_m'].items()))
    return '\n'.join(out)


def e2(model='smean2', fast_only=True, z_cut_km=300.):
    """Traction toward Earth's flow in the torque law on Bird's plates: the pass-45 table row for this density (PARAMS
    coefficients with traction, and a refit without suction / coupling). Prints force_balance_earth's report and the
    co-motion measures of trench_pairs."""
    from planet_reviews import force_balance_earth as fbe, trench_pairs as tp
    lmax, H, V = earth_flow(model, fast_only, z_cut_km)
    g, pid, names, cont = fbe.earth_plates(217); K = len(names)
    vfs = mf.gradient(V, lmax, g.xyz)                     # m/s
    st = fbe.earth_steps(names, fbe._rings()); T = fbe.terms(g, pid, cont, st, K); Am = fbe.assemble(T, K)
    met = fbe.Metric(g, pid, K)
    R = fbe.R
    tq = np.cross(g.xyz, vfs)*R*g.area[:, None]
    F0 = np.stack([tq[pid == k].sum(0) for k in range(K)]); F1 = np.stack([(tq*cont[:, None])[pid == k].sum(0) for k in range(K)])
    Am['F0'] = F0.reshape(-1); Am['F1'] = F1.reshape(-1)

    def solve(prm, trac=True):
        Cb, beta, Fsp, Cs, Fsu, Cc = prm
        A = Cb*(Am['D']+beta*Am['K'])+Cs*Am['Z']+Cc*Am['X']
        b = Fsp*Am['S']+fbe.F_RP*Am['P']+Fsu*Am['U']+(Cb*(Am['F0']+beta*Am['F1']) if trac else 0)
        w = np.linalg.solve(A, b)
        return {k: w[3*k:3*k+3] for k in range(K)}
    w_obs = fbe.earth_omega(g, pid, names)
    E = tp.earth_segments(); oo = E['io'] | E['ocb']; x = st['x']; L = E['L']; cls = E['cls']; ms = R*1e3*YR_S

    def pairs(w):
        va = np.cross(np.stack([w[int(k)] for k in st['A']]), x)*ms; vb = np.cross(np.stack([w[int(k)] for k in st['B']]), x)*ms
        close = -np.einsum('ij,ij->i', vb-va, st['nA']); rho = tp.rho(va, vb)
        row = [f"oo {tp.share(rho[oo] < .5, L[oo]):.2f} (<20: {tp.share(close[oo] < 20, L[oo]):.2f})"]+[
            f"{c} {tp.share(rho[cls == c] < .5, L[cls == c]):.2f}" for c in ('OSR', 'CRB', 'OTF', 'CTF', 'SUB', 'OCB', 'CCB')]
        wn = np.stack([w[int(k)] for k in pid]); coh = tp.coherence(g, pid, np.cross(wn, g.xyz), np.random.default_rng(0))
        print('    rho<0.5: '+'; '.join(row)+'\n    C(d): '+tp.coh_txt(coh))
    prm0 = fbe.params_coefficients()
    w = met.nnr(solve(prm0))
    fbe.report({'g': g, 'pid': pid, 'names': names, 'st': st, 'w': w, 'w_obs': w_obs, 'K': K, 'prm': prm0}, '  PARAMS + traction')
    pairs(w)
    w = met.nnr(solve(prm0, False))
    fbe.report({'g': g, 'pid': pid, 'names': names, 'st': st, 'w': w, 'w_obs': w_obs, 'K': K, 'prm': prm0}, '  PARAMS, traction off')
    pairs(w)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--traction', action='store_true')
    a = ap.parse_args(argv)
    print(e1_text(e1()))
    if a.traction:
        for key in VARIANTS[:3]:
            print(f'\nE2: {key}')
            e2(*key)


if __name__ == '__main__':
    main()
