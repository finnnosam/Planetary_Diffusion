"""Pass 52, reported: the circles against E3 split by degree band (l 1-3 and l 4-31), joint regression with a
10-degree-block bootstrap (seed 52, 2,000), and each band's std against the data's. Run from the repo root."""
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_tectonics import mantle_flow as mf
from planet_tectonics.params import PARAMS as PR
HD = DATA/'residual_topography'/'Residual depthtopographic measurements'/'Residual depthtopographic measurements'/'Holdt_et_al_2022_SupportingInformation'
PC = OUTPUTS/'planet-calibration'; A = 6371e3; RM, RW = PR['rho_mantle'], PR['rho_water']
xyz = lambda lo, la: np.stack([np.cos(np.radians(la))*np.cos(np.radians(lo)), np.cos(np.radians(la))*np.sin(np.radians(lo)), np.sin(np.radians(la))], 1)
E = np.load(PC/'precheck-p50-geoid'/'e3_loads.npz'); RAD = E['radii']; L = int(E['lmax']); z = (A-RAD)/1e3
base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08'); K = mf.FlowKernels(A, RAD, base, eta, L)
a = np.loadtxt(HD/'oceanic_measurements_circles_1_deg.xyz'); lo, la, d = a[:, 0], a[:, 1], a[:, 2]; X = xyz(lo, la)
blk = (np.floor((la+90)/10)*36+np.floor(np.mod(lo, 360)/10)).astype(int); ub, inv = np.unique(blk, return_inverse=True)
mem = [np.flatnonzero(inv == i) for i in range(ub.size)]; rng = np.random.default_rng(52)
boots = [np.concatenate([mem[i] for i in rng.integers(0, ub.size, ub.size)]) for _ in range(2000)]
li = mf.alm_l(L)
for name, f in (('thermo', mf.thermo_ratio(z)), ('alpha(z)', mf.alpha_ratio(z))):
    H = K.surface(E['loads_0']*f[:, None])[0]; H[li == 0] = 0
    e13 = mf.synthesis(np.where(li <= 3, H, 0), L, X)/(RM-RW)/1e3; e4 = mf.synthesis(np.where(li >= 4, H, 0), L, X)/(RM-RW)/1e3
    def fit(i):
        M = np.stack([np.ones(i.size), e13[i], e4[i]], 1); return np.linalg.lstsq(M, d[i], rcond=None)[0][1:]
    c = fit(np.arange(d.size)); bt = np.array([fit(i) for i in boots]); q = np.percentile(bt, [2.5, 97.5], 0)
    r13 = np.corrcoef(e13, d)[0, 1]; r4 = np.corrcoef(e4, d)[0, 1]
    print(f'{name}: std at circles e l1-3 {e13.std():.3f}, e l4-31 {e4.std():.3f}, d {d.std():.3f} km | joint slopes: l1-3 {c[0]:.2f} '
          f'[{q[0,0]:.2f}, {q[1,0]:.2f}], l4-31 {c[1]:.2f} [{q[0,1]:.2f}, {q[1,1]:.2f}] | r alone l1-3 {r13:+.3f}, l4-31 {r4:+.3f}')
