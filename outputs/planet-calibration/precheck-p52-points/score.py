"""Pass 52: E3 sampled at Holdt et al.'s (2022) oceanic residual-depth points (rule P52_RULE.txt). Run from the repo root.
argv: [nboot] [nnull]"""
import sys, json
import numpy as np, h5py
from planet_grid.paths import DATA, OUTPUTS
from planet_tectonics import mantle_flow as mf
from planet_tectonics.params import PARAMS as PR
NB = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
NN = int(sys.argv[2]) if len(sys.argv) > 2 else 500
HD = DATA/'residual_topography'/'Residual depthtopographic measurements'/'Residual depthtopographic measurements'/'Holdt_et_al_2022_SupportingInformation'
PC = OUTPUTS/'planet-calibration'
A = 6371e3; RHO_M, RHO_W = PR['rho_mantle'], PR['rho_water']

def xyz(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], 1)

def load(n):
    a = np.loadtxt(HD/f'oceanic_measurements_{n}_1_deg.xyz')
    return a[:, 0], a[:, 1], a[:, 2], -a[:, 3], a[:, 4]

E = np.load(PC/'precheck-p50-geoid'/'e3_loads.npz'); RAD = E['radii']; L = int(E['lmax']); z = (A-RAD)/1e3
base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08')
K = mf.FlowKernels(A, RAD, base, eta, L)
E15 = np.load(PC/'precheck-p51-mass'/'e3_area_15.npz')
assert np.allclose(E15['radii'], RAD)
models = {'thermo': K.surface(E['loads_0']*mf.thermo_ratio(z)[:, None])[0],
          'alpha(z)': K.surface(E['loads_0']*mf.alpha_ratio(z)[:, None])[0],
          'thermo, 15 mm/yr': K.surface(E15['loads_0']*mf.thermo_ratio(z)[:, None])[0]}
for H in models.values():
    H[mf.alm_l(L) == 0] = 0.

lon, lat, d, em, ep = load('circles'); X = xyz(lon, lat); err = 0.5*(em+ep)
blk = (np.floor((lat+90)/10)*36+np.floor(np.mod(lon, 360)/10)).astype(int)
ub, binv = np.unique(blk, return_inverse=True); members = [np.flatnonzero(binv == i) for i in range(ub.size)]
rng = np.random.default_rng(52)
boots = [np.concatenate([members[i] for i in rng.integers(0, ub.size, ub.size)]) for _ in range(NB)]
rots = []
for _ in range(NN):
    q, r = np.linalg.qr(rng.normal(size=(3, 3))); q = q*np.sign(np.diag(r)); rots.append(q if np.linalg.det(q) > 0 else -q)

def stats(e, d, w=None):
    w = np.ones_like(d) if w is None else w
    me, md = np.average(e, weights=w), np.average(d, weights=w)
    ve, vd = np.average((e-me)**2, weights=w), np.average((d-md)**2, weights=w)
    c = np.average((e-me)*(d-md), weights=w)
    return dict(P=np.sqrt(ve/vd), r=c/np.sqrt(ve*vd), s=c/ve, b=c/vd, me=me, md=md, sde=np.sqrt(ve), sdd=np.sqrt(vd))

# Holdt's own l 1-30 fit and the age grid
with h5py.File(HD/'degree_1_to_30.grd', 'r') as h:
    G30 = h['z'][:]; glon, glat = h['lon'][:], h['lat'][:]
g30 = G30[np.rint(lat-glat[0]).astype(int), np.rint(np.mod(lon, 360)-glon[0]).astype(int)]
with h5py.File(HD/'age.2020.1.GTS2012.2m_EDIT.grd', 'r') as h:
    alat, alon = h['lat'][:], h['lon'][:]
    ii = np.arange(15, alat.size, 30); jj = np.arange(15, alon.size, 30)          # 1-degree cell centres
    age1 = h['z'][ii[0]:ii[-1]+1:30, jj[0]:jj[-1]+1:30]
LA, LO = np.meshgrid(alat[ii], alon[jj], indexing='ij'); ocean = np.isfinite(age1)
XG = xyz(LO.ravel(), LA.ravel()); wg = np.cos(np.radians(LA.ravel()))

out = {}
for name, H in models.items():
    e = mf.synthesis(H, L, X)/(RHO_M-RHO_W)/1e3                     # km, water-loaded
    S = stats(e, d); Sw = stats(e, d, 1/err**2)
    bt = np.array([[stats(e[i], d[i])[k] for k in ('P', 'r', 's', 'b')] for i in boots])
    lo95, hi95 = np.percentile(bt, 2.5, 0), np.percentile(bt, 97.5, 0)
    nullr = np.array([stats(mf.synthesis(H, L, X@R.T)/(RHO_M-RHO_W)/1e3, d)['r'] for R in rots]) if NN else np.zeros(1)
    hg = mf.synthesis(H, L, XG)/RHO_M/1e3                                # km, air-loaded, on the 1-degree grid
    rms_glob = np.sqrt(np.average(hg**2, weights=wg)); rms_oc = np.sqrt(np.average(hg[ocean.ravel()]**2, weights=wg[ocean.ravel()]))
    rms_alm = np.sqrt(mf.spectrum(H/RHO_M, L)[1:].sum())/1e3
    Sg = stats(e, g30)
    tri = {}
    for n in ('down_triangles', 'up_triangles'):
        lo_, la_, v, m_, p_ = load(n); et = mf.synthesis(H, L, xyz(lo_, la_))/(RHO_M-RHO_W)/1e3
        res = {}
        for lab, mod in (('zero model', np.full_like(v, S['md'])), ('e mean-matched', et-S['me']+S['md']),
                         ('s*e', S['md']+S['s']*(et-S['me']))):
            ok = mod <= v+p_ if n == 'down_triangles' else mod >= v-m_
            res[lab] = float(ok.mean())
        tri[n] = res
    out[name] = dict(n=int(d.size), blocks=int(ub.size), **{k: float(v) for k, v in S.items()},
                     P_ci=[float(lo95[0]), float(hi95[0])], r_ci=[float(lo95[1]), float(hi95[1])],
                     s_ci=[float(lo95[2]), float(hi95[2])], b_ci=[float(lo95[3]), float(hi95[3])],
                     weighted={k: float(Sw[k]) for k in ('P', 'r', 's')},
                     null_r_p95=float(np.percentile(nullr, 95)), null_r_mean=float(nullr.mean()),
                     null_frac_ge=float((nullr >= S['r']).mean()),
                     rms_air_alm=float(rms_alm), rms_air_grid=float(rms_glob), rms_air_ocean=float(rms_oc),
                     vs_holdt_l30=dict(P=float(Sg['P']), r=float(Sg['r']), s=float(Sg['s']), sd_grid=float(Sg['sdd'])),
                     triangles=tri)
    o = out[name]
    print(f"{name}: n {o['n']} in {o['blocks']} blocks | P {o['P']:.2f} [{o['P_ci'][0]:.2f}, {o['P_ci'][1]:.2f}] "
          f"(std e {o['sde']:.3f} / d {o['sdd']:.3f} km; means e {o['me']:+.3f}, d {o['md']:+.3f})")
    print(f"   r {o['r']:+.3f} [{o['r_ci'][0]:+.2f}, {o['r_ci'][1]:+.2f}], null mean {o['null_r_mean']:+.3f} p95 {o['null_r_p95']:+.3f}, "
          f"null >= r {o['null_frac_ge']:.3f} | bracket b {o['b']:.2f} [{o['b_ci'][0]:.2f},{o['b_ci'][1]:.2f}] <= k <= 1/s {1/o['s']:.2f} "
          f"(s {o['s']:.2f} [{o['s_ci'][0]:.2f},{o['s_ci'][1]:.2f}])")
    print(f"   weighted P {o['weighted']['P']:.2f} r {o['weighted']['r']:+.3f} s {o['weighted']['s']:.2f} | air-loaded rms: alm {o['rms_air_alm']*1e3:.0f} m, "
          f"1-deg grid {o['rms_air_grid']*1e3:.0f} m, ocean {o['rms_air_ocean']*1e3:.0f} m | vs Holdt l1-30 at circles: P {Sg['P']:.2f} r {Sg['r']:+.3f} (grid std {Sg['sdd']:.3f})")
    print('   triangles met (zero model / e / s*e): ' + '; '.join(f"{n} " + ' / '.join(f'{x:.3f}' for x in r.values()) for n, r in tri.items()))
json.dump(out, open(PC/'precheck-p52-points'/'score.json', 'w'), indent=1)
