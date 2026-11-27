"""Pass 81 (journal 105): Earth's seismic lithosphere against its crust, by distance from the continent's edge.
On lithosphere.earth_continent's domain (CRUST1.0 continental types, total >= 20 km, ice sheets out; n = 217 nodes),
Hoggard et al. 2020's SL2013sv L. Per DIST_BINS_KM row: area-weighted L quantiles, and the correlation of the normal
scores of crust and L (a Gaussian copula's rho). Self-check: the copula join applied to Earth's own (distance, crust)
with an independent uniform rank v reproduces Earth's L by crust bin? Run from the repo root."""
import sys, numpy as np
from statistics import NormalDist
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_earth.reference import ice_mask
from planet_landscape import graph_landscape as gl
from planet_tectonics import lithosphere as lt
N = NormalDist(); inv = np.vectorize(N.inv_cdf); cdf = np.vectorize(N.cdf)
g = CubeSphere(217)
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
c1 = DATA/'crust1'
T = np.loadtxt(c1/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]
S = np.loadtxt(c1/'sedthk.xyz')[:, 2].reshape(180, 360)[r, c]
tn = np.array([l.split() for l in (c1/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
whole = ~np.isin(tn, lt.OCEANIC_TYPES) & (T >= 20)
d = gl.distance_to(~whole, g.nbr, g.dist)
cont = whole & ~ice_mask(g)
dd = np.loadtxt(DATA.joinpath(*lt.HOGGARD_SL2013SV), delimiter=',')
grid = np.full((361, 721), np.nan)
grid[np.round((90-dd[:, 1])*2).astype(int), np.round((dd[:, 0] % 360)*2).astype(int) % 720] = dd[:, 2]
grid[:, 720] = grid[:, 0]
L = grid[np.round((90-lat)*2).astype(int), np.round((lon % 360)*2).astype(int)]
m = cont & np.isfinite(L)
a, dk, Tt, Tc, LL = g.area[m], d[m]/1e3, T[m], T[m]-S[m], L[m]
print(f'Earth continent nodes {m.sum()} (of {cont.sum()}); seismic L area mean {np.average(LL, weights=a):.1f} km, '
      f'clipped to the lamb frame [93.78, 177.74] {np.average(np.clip(LL, 93.78, 177.74), weights=a):.1f}; frame ref {lt.LAMB_REF_KM}')
def nscore(x, w):
    o = np.argsort(x, kind='stable'); cw = np.cumsum(w[o]); u = np.empty(x.size); u[o] = (cw-.5*w[o])/cw[-1]
    return inv(np.clip(u, 1e-6, 1-1e-6))
Q = np.linspace(0, 1, 21)
B = lt.DIST_BINS_KM
for name, TX in (('total', Tt), ('crystalline', Tc)):
    print(f'\n== crust as {name}')
    rows, rhos = [], []
    for lo, hi in zip(B[:-1], B[1:]):
        k = (dk >= lo) & (dk < hi)
        zt, zl = nscore(TX[k], a[k]), nscore(LL[k], a[k])
        rho = float(np.sum(a[k]*zt*zl)/np.sqrt(np.sum(a[k]*zt**2)*np.sum(a[k]*zl**2)))
        q = lt._wquantile(LL[k], a[k], Q)
        rows.append(q); rhos.append(rho)
        print(f'  {lo:5.0f}-{hi:5.0f} km: nodes {k.sum():6d} area {a[k].sum()/a.sum():.3f} rho {rho:.3f} L p10/50/90 '
              f'{q[2]:.0f}/{q[10]:.0f}/{q[18]:.0f} crust p50 {np.median(TX[k]):.1f}')
    # self-check: join on Earth's own (distance, crust), independent v
    rng = np.random.default_rng(0)
    for trial, v in (('independent v', rng.uniform(size=dk.size)),):
        Lj = np.empty(dk.size)
        for i, (lo, hi) in enumerate(zip(B[:-1], B[1:])):
            k = (dk >= lo) & (dk < hi)
            z = rhos[i]*nscore(TX[k], a[k])+np.sqrt(1-rhos[i]**2)*inv(np.clip(v[k], 1e-6, 1-1e-6))
            Lj[k] = np.interp(cdf(z), Q, rows[i])
        print(f'  self-check ({trial}): L area mean by total-crust bin, Earth seismic / joined / joined clipped to frame')
        for lo, hi in [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]:
            k = (Tt >= lo) & (Tt < hi)
            print(f'    {lo}-{hi}: {np.average(LL[k], weights=a[k]):.0f} / {np.average(Lj[k], weights=a[k]):.0f} / '
                  f'{np.average(np.clip(Lj[k], 93.78, 177.74), weights=a[k]):.0f}')
        print(f'    dL/dT node least squares: Earth {np.polyfit(Tt, LL, 1)[0]:.2f}, joined {np.polyfit(Tt, Lj, 1)[0]:.2f} km/km')
    # by crust only (one row, no distance): rho
    zt, zl = nscore(TX, a), nscore(LL, a)
    print(f'  whole continent rho {np.sum(a*zt*zl)/np.sqrt(np.sum(a*zt**2)*np.sum(a*zl**2)):.3f}')
