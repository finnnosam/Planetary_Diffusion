"""D's flow kernels against an independent solution (forty-eighth pass; the reference owed since pass 45).

The reference shares no code or convention with planet_tectonics.mantle_flow (no propagator, no A matrix, no state
vector): the radial Stokes equations for one degree, derived from the spherical stress divergence (sympy, journal 69;
u_r = U P_l, u_th = V dP_l/dth, p = Pi P_l), solved in primitive variables by Chebyshev collocation, one domain per
constant-viscosity layer, split at the load:
  continuity  r U' + 2U - L V = 0
  r-momentum  2 eta U'' + 4 eta U'/r - (L+4) eta U/r^2 - L eta V'/r + 3 L eta V/r^2 - Pi' = 0
  th-momentum eta V'' + 2 eta V'/r - 2 L eta V/r^2 + eta U'/r + 4 eta U/r^2 - Pi/r = 0
with U, V, sigma_rr = -Pi + 2 eta U', sigma_rt = eta (r V' + U - V)/r continuous between layers, sigma_rr jumping by
g sigma at the load (the body force -g drho integrated across the sheet), free slip at both boundaries. The matrix is
equilibrated (thin layers and viscosity ratios to ~6,000 make it badly scaled otherwise)."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_tectonics import mantle_flow as mf

A_M = 6371e3
VISC = DATA/'mantle'/'visc.sh08'


def _cheb(N):
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2, np.ones(N-1), 2])*(-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    D = np.outer(c, 1/c)/(X-X.T+np.eye(N+1))
    return D-np.diag(D.sum(1)), x


def reference_kernel(l, load_r, edges, etas, N=16, g=mf.G, a=A_M, eta0=1e21):
    """(KH, KV) as FlowKernels defines them: -sigma_rr(a)/g and u_h(a) per unit sheet (1 kg/m^2) at non-dimensional
    radius load_r; edges: layer boundaries (non-dimensional radius, CMB to 1), etas: their viscosities (Pa s).
    Units: r in a, eta in eta0, stress in g sigma (so the jump is 1), velocity in g a / eta0."""
    edges, etas = list(edges), list(etas)
    k = int(np.searchsorted(edges, load_r))-1
    if min(abs(edges[k]-load_r), abs(edges[k+1]-load_r)) > 1e-12:
        edges.insert(k+1, load_r); etas.insert(k, etas[k]); split = k+1
    else:
        split = int(np.argmin(np.abs(np.array(edges)-load_r)))
    nd, n, L = len(edges)-1, N+1, l*(l+1.)
    D0, x0 = _cheb(N)
    nu = 3*n*nd
    A = np.zeros((nu, nu)); b = np.zeros(nu)
    sl = lambda d, q: slice(3*n*d+q*n, 3*n*d+(q+1)*n)          # q: 0 U, 1 V, 2 Pi
    ops, row = [], 0
    for d in range(nd):
        r0, r1 = edges[d], edges[d+1]
        r = r0+(x0[::-1]+1)/2*(r1-r0)
        Dm = D0[::-1, ::-1]*2/(r1-r0); D2 = Dm@Dm
        e = etas[d]/eta0
        I = np.eye(n); R = np.diag(r); Ri = np.diag(1/r); Ri2 = np.diag(1/r**2)
        blocks = [((R@Dm+2*I, -L*I, 0*I), range(n)),
                  ((2*e*D2+4*e*Ri@Dm-(L+4)*e*Ri2, -L*e*Ri@Dm+3*L*e*Ri2, -Dm), range(1, n-1)),
                  ((e*Ri@Dm+4*e*Ri2, e*D2+2*e*Ri@Dm-2*L*e*Ri2, -Ri), range(1, n-1))]
        for (bu, bv, bp), js in blocks:
            for j in js:
                A[row, sl(d, 0)] = bu[j]; A[row, sl(d, 1)] = bv[j]; A[row, sl(d, 2)] = bp[j]; row += 1
        ops.append({'srr': (2*e*Dm, -I), 'srt': (e*Ri, e*(Dm-Ri))})

    def srr(d, j):
        v = np.zeros(nu); v[sl(d, 0)] = ops[d]['srr'][0][j]; v[sl(d, 2)] = ops[d]['srr'][1][j]; return v

    def srt(d, j):
        v = np.zeros(nu); v[sl(d, 0)] = ops[d]['srt'][0][j]; v[sl(d, 1)] = ops[d]['srt'][1][j]; return v

    def unit(d, q, j):
        v = np.zeros(nu); v[sl(d, q).start+j] = 1.; return v
    for v in (unit(0, 0, 0), srt(0, 0), unit(nd-1, 0, n-1), srt(nd-1, n-1)):
        A[row] = v; row += 1
    for d in range(nd-1):
        for v in (unit(d, 0, n-1)-unit(d+1, 0, 0), unit(d, 1, n-1)-unit(d+1, 1, 0), srt(d, n-1)-srt(d+1, 0)):
            A[row] = v; row += 1
        A[row] = srr(d+1, 0)-srr(d, n-1); b[row] = 1. if d+1 == split else 0.; row += 1
    assert row == nu
    dr = np.ones(nu); dc = np.ones(nu)
    for _ in range(30):
        rr = np.sqrt(np.abs(A).max(1)); cc = np.sqrt(np.abs(A).max(0)); rr[rr == 0] = 1; cc[cc == 0] = 1
        A = A/rr[:, None]/cc[None, :]; dr *= rr; dc *= cc
    x = np.linalg.solve(A, b/dr)/dc
    return -srr(nd-1, n-1)@x, x[sl(nd-1, 1)][n-1]*g*a/eta0


DEPTHS_KM = (325., 575., 1075., 1525., 1925., 2525., 2825.)   # 1925 and 1075 lie 5 km from sh08 layer boundaries
DEGREES = (1, 2, 3, 8, 20)


def _check(base, eta, edges):
    radii = A_M-np.array(DEPTHS_KM)*1e3
    K = mf.FlowKernels(A_M, radii, base, eta, max(DEGREES))
    for l in DEGREES:
        vmax = np.abs(K.KV[l]).max()
        for j, r in enumerate(radii):
            kh, kv = reference_kernel(l, r/A_M, edges, eta)
            assert abs(kh-K.KH[l, j]) < 2e-4, (l, DEPTHS_KM[j], kh, K.KH[l, j])
            if abs(K.KV[l, j]) > 1e-3*vmax:
                assert abs(kv/K.KV[l, j]-1) < 2e-3, (l, DEPTHS_KM[j], kv, K.KV[l, j])


def test_kernels_match_the_reference_isoviscous():
    c = mf.CMB_RADIUS_M/A_M
    _check(np.array([c]), np.array([1e21]), [c, 1.])


@pytest.mark.skipif(not VISC.exists(), reason='data/mantle/visc.sh08 not present')
def test_kernels_match_the_reference_sh08():
    base, eta = mf.read_viscosity(VISC)
    _check(base, eta, list(np.append(np.append(mf.CMB_RADIUS_M/A_M, base[1:]), 1.)))
