"""Pass 50: D's CMB kernel KC and geoid kernel against the independent collocation reference (tests/test_flow_reference),
every degree 1-31 at the 51 load depths, for sh08 and an isoviscous mantle; plus limits."""
import sys, numpy as np
sys.path.insert(0, 'tests')
from test_flow_reference import reference_kernel
from planet_tectonics import mantle_flow as mf
from planet_grid.paths import DATA
a = 6371e3; c = mf.CMB_RADIUS_M
radii = mf.load_radii(a, 300e3)
for name, (base, eta) in (('sh08', mf.read_viscosity(DATA/'mantle'/'visc.sh08')), ('isoviscous', (np.array([c/a]), np.array([1e21])))):
    edges = list(np.append(np.append(c/a, base[1:]), 1.))
    K = mf.FlowKernels(a, radii, base, eta, 31)
    KN, _ = K.geoid_kernels()
    wc = wh = wn = 0.
    for l in range(1, 32):
        kk = 4*np.pi*mf.GRAV_CONST/((2*l+1)*mf.G)
        for j, r in enumerate(radii):
            kh, kv, kc = reference_kernel(l, r/a, edges, eta)
            wc = max(wc, abs(kc-K.KC[l, j])); wh = max(wh, abs(kh-K.KH[l, j]))
            nref = kk*a*((r/a)**(l+2)+kh+kc*(c/a)**(l+2))
            if l > 1:
                wn = max(wn, abs(nref-KN[l, j])/np.abs(KN[l]).max())
    print(f'{name}: max |KC_D - KC_ref| {wc:.1e}; max |KH_D - KH_ref| {wh:.1e}; max |KN_D - KN_ref| / degree max {wn:.1e}')
    for l in (2, 3, 4, 8):
        print(f'  l={l}: geoid kernel (m per kg/m^2) at z = ' + ', '.join(f'{(a-radii[j])/1e3:.0f} km {KN[l, j]:+.2e}' for j in (0, 7, 14, 24, 34, 44, 50))
              + f';  KC at 2825 km {K.KC[l, -1]:+.3f}, KH+KC(c/a)^(l+2)+(r/a)^(l+2) at 1025 km {KN[l, 14]/(4*np.pi*mf.GRAV_CONST/((2*l+1)*mf.G)*a):+.3f}')
    # self-gravitation off by zero densities equals the plain kernel
    KN0, KH0 = K.geoid_kernels(True, 0., 0.)
    print('  SG with zero densities = plain:', float(np.abs(KN0-KN)[2:].max()), float(np.abs(KH0-K.KH)[1:].max()))
# limits: sheets at the boundaries are compensated (no geoid)
base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08')
K = mf.FlowKernels(a, np.array([a*(1-1e-7), c*(1+1e-7)]), base, eta, 8)
KN, _ = K.geoid_kernels()
print('limits: surface sheet KH', K.KH[2:5, 0].round(4), 'KC', K.KC[2:5, 0].round(4), '| CMB sheet KH', K.KH[2:5, 1].round(4), 'KC', K.KC[2:5, 1].round(4),
      '| N (m per kg/m^2)', KN[2:5].ravel())
