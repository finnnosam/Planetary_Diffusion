"""Pass 52, reported diagnostics after score.py (not in the rule): Holdt's l 1-30 fit against E3 over the ocean and by
degree, the circles by ocean basin and by E3's value. Run from the repo root."""
import numpy as np, h5py
from planet_grid.paths import DATA, OUTPUTS
from planet_tectonics import mantle_flow as mf
from planet_tectonics.params import PARAMS as PR
HD = DATA/'residual_topography'/'Residual depthtopographic measurements'/'Residual depthtopographic measurements'/'Holdt_et_al_2022_SupportingInformation'
PC = OUTPUTS/'planet-calibration'; A = 6371e3; RM, RW = PR['rho_mantle'], PR['rho_water']
xyz = lambda lo, la: np.stack([np.cos(np.radians(la))*np.cos(np.radians(lo)), np.cos(np.radians(la))*np.sin(np.radians(lo)), np.sin(np.radians(la))], 1)
E = np.load(PC/'precheck-p50-geoid'/'e3_loads.npz'); RAD = E['radii']; L = int(E['lmax']); z = (A-RAD)/1e3
base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08'); K = mf.FlowKernels(A, RAD, base, eta, L)
H = K.surface(E['loads_0']*mf.thermo_ratio(z)[:, None])[0]; H[mf.alm_l(L) == 0] = 0
with h5py.File(HD/'degree_1_to_30.grd', 'r') as h:
    G = h['z'][:].astype(float); glon, glat = h['lon'][:], h['lat'][:]
G = G[:, :-1]; LA, LO = np.meshgrid(glat, glon[:-1], indexing='ij')           # drop the repeated 360 column
with h5py.File(HD/'age.2020.1.GTS2012.2m_EDIT.grd', 'r') as h:
    alat, alon = h['lat'][:], h['lon'][:]
    ia = np.rint((LA.ravel()-alat[0])*30).astype(int); ja = np.rint((((LO.ravel()+180) % 360)-180-alon[0])*30).astype(int)
    Z = h['z']; age = np.array([Z[i, j] for i, j in zip(ia, ja)]) if False else None
    rows = np.unique(ia); age = np.full(ia.size, np.nan)
    for i in rows:
        m = ia == i; age[m] = Z[i, :][ja[m]]
oc = np.isfinite(age); w = np.cos(np.radians(LA.ravel())); X = xyz(LO.ravel(), LA.ravel()); g = G.ravel()
ew = mf.synthesis(H, L, X)/(RM-RW)/1e3; ea = mf.synthesis(H, L, X)/RM/1e3
rms = lambda v, m: np.sqrt(np.average(v[m]**2, weights=w[m]))
cor = lambda a, b, m: np.average(a[m]*b[m], weights=w[m])/rms(a, m)/rms(b, m)
print(f'ocean fraction (age defined) {np.average(oc, weights=w):.3f}')
print(f'Holdt l1-30 rms: global {rms(g, oc|~oc):.3f} km, ocean {rms(g, oc):.3f}, land {rms(g, ~oc):.3f}')
print(f'E3 thermo rms: ocean water-loaded {rms(ew, oc):.3f} km (air {rms(ea, oc):.3f}), land air {rms(ea, ~oc):.3f}')
print(f'ocean: E3 water / Holdt {rms(ew, oc)/rms(g, oc):.2f}, corr {cor(ew, g, oc):+.3f}; land: E3 air / Holdt {rms(ea, ~oc)/rms(g, ~oc):.2f}, corr {cor(ea, g, ~oc):+.3f}; global corr {cor(ea, g, oc|~oc):+.3f}')
# degree split by quadrature (1-degree grid, cos-lat weights)
dA = (np.pi/180)**2*w
for name, f in (('Holdt l1-30', g), ('E3 air-loaded', ea)):
    alm = mf.point_alm(X, f*dA, L); P = mf.spectrum(alm, L)
    print(f'{name}: rms l1-31 {np.sqrt(P[1:].sum()):.3f} km, share l1-3 {P[1:4].sum()/P[1:].sum():.2f}, l4-8 {P[4:9].sum()/P[1:].sum():.2f}, l9-31 {P[9:].sum()/P[1:].sum():.2f}; rms l1-3 {np.sqrt(P[1:4].sum()):.3f}')
# circles by basin and by E3 value
a = np.loadtxt(HD/'oceanic_measurements_circles_1_deg.xyz'); lo, la, d = a[:, 0], a[:, 1], a[:, 2]
e = mf.synthesis(H, L, xyz(lo, la))/(RM-RW)/1e3; l360 = np.mod(lo, 360)
basins = {'Pacific (120-290 E)': (l360 >= 120) & (l360 < 290) & (la > -60),
          'Atlantic (290-20 E)': ((l360 >= 290) | (l360 < 20)) & (la > -60),
          'Indian (20-120 E)': (l360 >= 20) & (l360 < 120) & (la > -60), 'Southern (< 60 S)': la <= -60}
for k, m in basins.items():
    if m.sum() > 5:
        c = np.cov(e[m], d[m]); print(f'{k:22s} n {m.sum():4d}: std e {e[m].std():.3f} d {d[m].std():.3f} P {e[m].std()/d[m].std():.2f} r {c[0,1]/np.sqrt(c[0,0]*c[1,1]):+.2f}; mean e {e[m].mean():+.3f} d {d[m].mean():+.3f}')
for lo_, hi_ in ((-9, -0.5), (-0.5, 0.), (0., 0.5), (0.5, 9)):
    m = (e >= lo_) & (e < hi_)
    print(f'E3 in [{lo_:+.1f}, {hi_:+.1f}) km: n {m.sum():4d}, mean e {e[m].mean():+.3f}, mean d {d[m].mean():+.3f}')
