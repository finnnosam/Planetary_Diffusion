"""Independent reference for D's flow kernels (pass 48): primitive-variable Chebyshev collocation of the radial Stokes
equations derived in derive.py (sympy, from the spherical stress divergence). Shares no code or convention with
planet_tectonics.mantle_flow: no propagator, no A matrix, no state vector.

Per degree l, unknowns U (radial velocity), V (u_th = V dP/dth), Pi (pressure) on Gauss-Lobatto nodes of each domain
(a constant-viscosity layer, split at the load radius). Equations (eta constant inside a domain, L = l(l+1)):
  continuity  r U' + 2U - L V = 0                                   (all nodes)
  r-momentum  2 eta U'' + 4 eta U'/r - (L+4) eta U / r^2 - L eta V'/r + 3 L eta V / r^2 - Pi' = 0    (interior)
  th-momentum eta V'' + 2 eta V'/r - 2 L eta V / r^2 + eta U'/r + 4 eta U / r^2 - Pi/r = 0            (interior)
Interfaces: U, V, sigma_rr = -Pi + 2 eta U', sigma_rt = eta (r V' + U - V)/r continuous; at the load radius
sigma_rr(+) - sigma_rr(-) = g sigma (the integral of the body force -g drho across the sheet).
Free slip at both boundaries: U = 0, sigma_rt = 0.
Non-dimensional: r in units of a, eta in units of eta0, stress in units of g sigma (so the jump is 1).
Returns KH = -sigma_rr(a)/g per unit sigma (kg/m^2 per kg/m^2) and KV = u_h(a) (m/s per kg/m^2), as D's FlowKernels.
"""
import numpy as np


def cheb(N):
    """Trefethen's Chebyshev differentiation matrix on x_j = cos(pi j / N), j = 0..N (x from 1 to -1)."""
    x = np.cos(np.pi*np.arange(N+1)/N)
    c = np.hstack([2, np.ones(N-1), 2])*(-1)**np.arange(N+1)
    X = np.tile(x, (N+1, 1)).T
    dX = X-X.T
    D = np.outer(c, 1/c)/(dX+np.eye(N+1))
    D = D-np.diag(D.sum(1))
    return D, x


def kernel(l, load_r, edges, etas, N=24, g=9.81, a=6371e3, eta0=1e21, return_profile=False):
    """edges: domain boundaries (non-dim radius, increasing, from the CMB to 1); etas: viscosity per layer (Pa s).
    load_r: non-dim load radius. Returns (KH, KV)."""
    edges = list(edges); etas = list(etas)
    # split the layer that holds the load
    k = np.searchsorted(edges, load_r)-1
    split = None
    if 0 < load_r < 1 and abs(edges[k]-load_r) > 1e-12 and abs(edges[k+1]-load_r) > 1e-12:
        edges.insert(k+1, load_r); etas.insert(k, etas[k]); split = k+1
    else:
        split = int(np.argmin(np.abs(np.array(edges)-load_r)))
    nd = len(edges)-1
    L = l*(l+1.)
    D0, x0 = cheb(N)
    n = N+1
    nu = 3*n*nd
    A = np.zeros((nu, nu)); b = np.zeros(nu)
    iu = lambda d: slice(3*n*d, 3*n*d+n)
    iv = lambda d: slice(3*n*d+n, 3*n*d+2*n)
    ip = lambda d: slice(3*n*d+2*n, 3*n*d+3*n)
    row = 0
    ops = []
    for d in range(nd):
        r0, r1 = edges[d], edges[d+1]
        r = r0+(x0[::-1]+1)/2*(r1-r0)                  # increasing radius
        Dm = D0[::-1, ::-1]*2/(r1-r0)                   # d/dr on the increasing grid
        D2 = Dm@Dm
        e = etas[d]/eta0
        R = np.diag(r); Ri = np.diag(1/r); Ri2 = np.diag(1/r**2); I = np.eye(n)
        # continuity at all nodes
        C_u = R@Dm+2*I; C_v = -L*I
        # momentum
        Mr_u = 2*e*D2+4*e*Ri@Dm-(L+4)*e*Ri2; Mr_v = -L*e*Ri@Dm+3*L*e*Ri2; Mr_p = -Dm
        Mt_u = e*Ri@Dm+4*e*Ri2; Mt_v = e*D2+2*e*Ri@Dm-2*L*e*Ri2; Mt_p = -Ri
        # traction operators (rows evaluated at the ends)
        Srr_u = 2*e*Dm; Srr_p = -I
        Srt_u = e*Ri; Srt_v = e*(Dm-Ri)
        ops.append(dict(r=r, Dm=Dm, e=e, Srr_u=Srr_u, Srr_p=Srr_p, Srt_u=Srt_u, Srt_v=Srt_v))
        for j in range(n):
            A[row, iu(d)] = C_u[j]; A[row, iv(d)] = C_v[j]; row += 1
        for j in range(1, n-1):
            A[row, iu(d)] = Mr_u[j]; A[row, iv(d)] = Mr_v[j]; A[row, ip(d)] = Mr_p[j]; row += 1
            A[row, iu(d)] = Mt_u[j]; A[row, iv(d)] = Mt_v[j]; A[row, ip(d)] = Mt_p[j]; row += 1

    def srr(d, j, s=1.):
        o = ops[d]
        v = np.zeros(nu); v[iu(d)] = s*o['Srr_u'][j]; v[ip(d)] = s*o['Srr_p'][j]; return v

    def srt(d, j, s=1.):
        o = ops[d]
        v = np.zeros(nu); v[iu(d)] = s*o['Srt_u'][j]; v[iv(d)] = s*o['Srt_v'][j]; return v

    def unit(sl, j, s=1.):
        v = np.zeros(nu); v[sl.start+j] = s; return v
    # boundaries: free slip
    A[row] = unit(iu(0), 0); row += 1
    A[row] = srt(0, 0); row += 1
    A[row] = unit(iu(nd-1), n-1); row += 1
    A[row] = srt(nd-1, n-1); row += 1
    # interfaces
    for d in range(nd-1):
        A[row] = unit(iu(d), n-1)-unit(iu(d+1), 0); row += 1
        A[row] = unit(iv(d), n-1)-unit(iv(d+1), 0); row += 1
        A[row] = srt(d, n-1)-srt(d+1, 0); row += 1
        A[row] = srr(d+1, 0)-srr(d, n-1)                 # sigma_rr(+) - sigma_rr(-) = jump
        if d+1 == split:
            b[row] = 1.
        row += 1
    if load_r >= 1-1e-12:                                # a sheet at the surface: the surface stress carries it
        pass
    assert row == nu, (row, nu)
    # equilibrate (Ruiz: alternate row / column max-norm scaling): thin domains and viscosity ratios to ~6,000 make
    # the raw collocation matrix badly scaled
    dr = np.ones(nu); dc = np.ones(nu); As = A.copy()
    for _ in range(30):
        rr = np.sqrt(np.abs(As).max(1)); cc = np.sqrt(np.abs(As).max(0))
        rr[rr == 0] = 1; cc[cc == 0] = 1
        As = As/rr[:, None]/cc[None, :]; dr = dr*rr; dc = dc*cc
    sol = np.linalg.solve(As, b/dr)/dc
    s_top = srr(nd-1, n-1)@sol                           # sigma_rr(a) in units of g sigma
    if load_r >= 1-1e-12:
        s_top = s_top-1.
    KH = -s_top                                          # -sigma_rr/g per unit sigma
    u_scale = g*a/eta0                                   # velocity unit for sigma = 1 kg/m^2
    KV = sol[iv(nd-1)][n-1]*u_scale
    if return_profile:
        return KH, KV, ops, sol
    return KH, KV
