"""Pass 80 (journal 104): Earth's continental lithosphere thickness against tectonothermal age, in Lamb's frame.

Hoggard et al. 2020 (Nat. Geosci. 13, 504; https://osf.io/twksd) LAB depth maps (0.5 deg) sampled at the n = 217
graph's nodes on Earth's continent as C defines it (lithosphere.earth_continent: CRUST1.0 continental types, total
>= 20 km, ice sheets out; ages CLASS_AGE_MA). Fits the conductive plate-cooling form

    L(t) = L_inf - (L_inf - L_a) exp(-t / tau)

by area-weighted least squares (tau on a grid, L_inf and L_a linear at each tau); the mechanism's own tau for a plate
of thickness L_inf is L_inf^2 / (pi^2 kappa), kappa 1e-6 m^2/s. Means, not medians: the isostatic term is linear in L.
Prints the law, the class means against it and against C's Artemieva line, the reference (area mean on Earth's ages),
L at 3.6 Ga, Earth's dL/dT on nodes, and the same fit on the three comparison models (LAB > 350 km dropped: the
3D2015-07Sv / SLNAAFSA outliers). Run from the repo root."""
import sys, numpy as np
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_earth.reference import ice_mask
from planet_landscape import graph_landscape as gl
from planet_tectonics.lithosphere import CLASS_AGE_MA, OCEANIC_TYPES

ROOT = DATA/'Hoggard2020'/'twksd-osfstorage-archive'/'supplementary_datasets'/'lithospheric_thickness_maps'/'global'
MODELS = ('SL2013sv', 'CAM2016', '3D2015-07Sv', 'SLNAAFSA')
g = CubeSphere(217)
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
c1 = DATA/'crust1'
T = np.loadtxt(c1/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]
tn = np.array([l.split() for l in (c1/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
whole = ~np.isin(tn, OCEANIC_TYPES) & (T >= 20)
cont = whole & ~ice_mask(g)
age = np.array([CLASS_AGE_MA.get(t, np.nan) for t in tn]); a = g.area


def sample(model):
    d = np.loadtxt(ROOT/model/f'{model}.csv', delimiter=',')
    grid = np.full((361, 721), np.nan)
    grid[np.round((90-d[:, 1])*2).astype(int), np.round((d[:, 0] % 360)*2).astype(int) % 720] = d[:, 2]
    grid[:, 720] = grid[:, 0]
    return grid[np.round((90-lat)*2).astype(int), np.round((lon % 360)*2).astype(int)]


def fit(t, L, w):
    best = None
    for tau in np.arange(10., 1500.1, 1.):
        A = np.vstack([1-np.exp(-t/tau), np.exp(-t/tau)]).T
        sw = np.sqrt(w); x = np.linalg.lstsq(A*sw[:, None], L*sw, rcond=None)[0]
        e = float((w*(A@x-L)**2).sum()/w.sum())
        if best is None or e < best[0]:
            best = (e, tau, x[0], x[1])
    return best


def law(t, Linf, La, tau):
    return Linf-(Linf-La)*np.exp(-np.asarray(t)/tau)


out = {}
for model in MODELS:
    L = sample(model)
    m = cont & np.isfinite(L) & (L <= 350)
    e, tau, Linf, La = fit(age[m], L[m], a[m])
    ref = float((law(age[m], Linf, La, tau)*a[m]).sum()/a[m].sum())
    mean_obs = float((L[m]*a[m]).sum()/a[m].sum())
    kappa_tau = (Linf*1e3)**2/(np.pi**2*1e-6)/3.15576e13
    w = a[m]; tm = T[m]; Lm = L[m]
    sl = np.polyfit(tm, Lm, 1, w=np.sqrt(w))[0]
    out[model] = (Linf, La, tau)
    print(f'== {model}: {m.sum()} continental nodes (of {cont.sum()}; dropped {(cont & ~m).sum()})')
    print(f'   law L = {Linf:.2f} - ({Linf:.2f} - {La:.2f}) exp(-t / {tau:.0f} Myr); rms {np.sqrt(e):.1f} km;'
          f' plate-cooling tau for L_inf: {kappa_tau:.0f} Myr')
    print(f'   area mean L observed {mean_obs:.2f} km, law on Earth ages (ref) {ref:.2f}; L(0) {La:.1f}, '
          f'L(3600) {law(3600, Linf, La, tau):.2f}; dL/dT on nodes {sl:.2f} km per km')
    if model == 'SL2013sv':
        print('   age (Ma) | area share | observed mean / p50 | law | Artemieva 93.6 + 0.04 t')
        for t0 in sorted(set(age[m])):
            k = m & (age == t0)
            print(f'   {t0:6.0f} | {a[k].sum()/a[m].sum():.3f} | {(L[k]*a[k]).sum()/a[k].sum():5.1f} / '
                  f'{np.median(L[k]):5.1f} | {law(t0, Linf, La, tau):5.1f} | {93.6+.04*t0:5.1f}')
        print('   crust T bin (km) | observed mean L | law mean L | area')
        for lo, hi in [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]:
            k = m & (T >= lo) & (T < hi)
            print(f'   {lo}-{hi} | {(L[k]*a[k]).sum()/a[k].sum():5.1f} | '
                  f'{(law(age[k], Linf, La, tau)*a[k]).sum()/a[k].sum():5.1f} | {a[k].sum()/a[m].sum():.3f}')
        k1 = m & (age <= 300); k2 = m & (age >= 1500)
        print(f'   contrast young (<= 300 Ma) vs old (>= 1.5 Ga) law means: '
              f'{(law(age[k2],Linf,La,tau)*a[k2]).sum()/a[k2].sum()-(law(age[k1],Linf,La,tau)*a[k1]).sum()/a[k1].sum():.1f} km'
              f' (Artemieva {(.04*age[k2]*a[k2]).sum()/a[k2].sum()-(.04*age[k1]*a[k1]).sum()/a[k1].sum():.1f})')
print('spread over models: L_inf %.1f-%.1f, L_a %.1f-%.1f, tau %.0f-%.0f' % (
    min(v[0] for v in out.values()), max(v[0] for v in out.values()), min(v[1] for v in out.values()),
    max(v[1] for v in out.values()), min(v[2] for v in out.values()), max(v[2] for v in out.values())))
print("Lamb's k 10.753 m per km x rho_mantle 3300 -> drho %.2f kg/m^3" % (10.753*3.3))
