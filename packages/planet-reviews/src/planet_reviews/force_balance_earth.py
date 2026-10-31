"""The torque-balance force law on Earth's plates (B, forty-third pass): driving forces and resistances from Bird's (2003)
boundaries, fitted to Earth's measured plate motions; the Earth test for planet_tectonics' `force_law` 'torque'.

  python -m planet_reviews.force_balance_earth [--n 217] [--json out.json]

Each plate's angular velocity w_i (rad/s) solves the balance of torques about Earth's centre (Forsyth & Uyeda 1975,
Geophys. J. R. astr. Soc. 43; no inertia):

  F_sp S_i + F_rp P_i + F_su U_i = [C_b (D_i + beta K_i) + C_s Z_i] w_i + C_c sum_j X_ij (w_i - w_j)

  S_i  slab pull: per unit force, sum over the trench steps where i is the downgoing plate of L R (x x n), n the
       boundary normal pointing out of i (Bird's slash names the subducting plate: A/B, B subducts; A\\B, A subducts)
  P_i  ridge push (OSR steps, both plates): sum L R (x x -n), the push away from the ridge
  U_i  slab suction (trench steps where i overrides): sum L R (x x n), toward the trench (Conrad & Lithgow-Bertelloni
       2002, Science 298); fitted, and dropped when it does not improve the fit
  D_i  basal drag under the plate: R^2 sum_nodes a (I - x x^T) (traction -C_b v, a Couette asthenosphere)
  K_i  the same over its continental nodes (keels; Earth: area above -1000 m, as earth-plate-speeds-nnr.json)
  Z_i  slab resistance: sum over i's downgoing trench steps of L R^2 q q^T, q = x x n (drag on the slab, linear in
       the plate's velocity toward its trench; Forsyth & Uyeda's 'slab resistance', which with slab pull sets a
       terminal speed independent of plate size)
  X_ij convergence without a slab (CCB and OCB steps): L R^2 q q^T on the pair's relative motion

F_rp is fixed at 3e12 N/m (half-space cooling of ~80-100 Myr lithosphere; Turcotte & Schubert 2014, Geodynamics
eq. 6.405; Parsons & Richter 1980): it sets the absolute scale. The other coefficients are fitted to Earth's motions
(Bird's PB2002 poles, NNR frame) by minimising the area-weighted velocity misfit: Earth measurement, not the planet.

Test bands (fixed before the first fit, forty-third pass):
  T-a  variance reduction of the NNR velocity field (1 - int |v_law - v_earth|^2 / int |v_earth|^2) >= 0.6
  T-b  slab-less large plates AF, AN, NA, EU, SA: law speed <= 3.5 cm/yr each (Earth 1.2-2.9)
  T-c  slab plates PA, NZ, CO, AU, PS: law / Earth speed within 1/1.35 ... 1.35
  T-d  boundary rates from the law's poles on Bird's steps (length-weighted p50): OSR opening and SUB closing within
       25% of Bird's; OCB closing <= 20 mm/yr; CCB closing within a factor 2
  T-e  plausibility: the fit is held inside BOX (C_b 1e13-3e15 Pa s/m: asthenosphere 1e18-1e21 Pa s over 100-300 km;
       F_sp 0.5e13-5e13 N/m: Turcotte & Schubert, Schellart 2004; terminal slab speed 2-30 cm/yr; keels x1.1-x33);
       T-e fails when the optimum sits on the box's edge in C_b or F_sp
The current length law (planet_tectonics' 'lengths', speed from lengths, pole from the normalised pull) is scored on
T-a ... T-d beside it.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_grid.jit import njit
from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import DATA, OUTPUTS
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_tectonics.params import PARAMS

R = 6.371e6
YR = 3.15576e7
F_RP = 3e12                                   # N/m, ridge push (fixed: sets the scale)

# Bird (2003) PB2002_poles.xls: plate, pole latitude, pole longitude (deg), rate (deg/Myr), Pacific frame
POLES = (('AF', 59.16, -73.174, 0.9270359), ('AM', 57.645, -83.736, 0.9309), ('AN', 64.315, -83.984, 0.8695),
         ('AP', 33.639, -81.177, 0.916), ('AR', 59.658, -33.193, 1.1616), ('AS', 74.275, -87.237, 0.6497),
         ('AT', 56.283, 8.932, 1.64), ('AU', 60.08, 1.742, 1.07438632), ('BH', 12.559, 87.957, 0.3029),
         ('BR', 45.9, -111.0, 0.2), ('BS', 16.007, 122.442, 2.125), ('BU', 8.894, -75.511, 2.667),
         ('CA', 54.313, -79.431, 0.904), ('CL', 10.13, -45.57, 0.309), ('CO', 36.823, -108.629, 1.9975018),
         ('CR', -12.628, 175.127, 3.605), ('EA', 28.3, 66.4, 11.4), ('EU', 61.066, -85.819, 0.8591457),
         ('FT', -10.158, -178.305, 4.848), ('GP', 9.399, 79.69, 5.275), ('IN', 60.494, -30.403, 1.10335918),
         ('JF', 35.0, 26.0, 0.506786), ('JZ', 35.91, 70.166, 22.52), ('KE', 47.521, -3.115, 2.831),
         ('MA', 43.777, 149.205, 1.278), ('MN', -3.037, 150.456, 51.3), ('MO', 59.589, 78.88, 0.8927),
         ('MS', 11.103, -56.746, 4.07), ('NA', 48.709, -78.167, 0.74860898), ('NB', -4.0, 139.0, 0.33),
         ('ND', 58.664, -89.003, 0.7009), ('NH', 13.0, -12.0, 2.7), ('NI', 6.868, -168.868, 3.255),
         ('NZ', 55.578, -90.096, 1.35990764), ('OK', 55.421, -82.859, 0.845), ('ON', 48.351, 142.415, 2.853),
         ('PA', 0.0, 0.0, 0.0), ('PM', 54.058, -90.247, 0.9069), ('PS', -1.2, -45.8, 1.0),
         ('RI', 26.7, -105.2, 4.6923076923), ('SA', 54.999, -85.752, 0.63654234), ('SB', 10.61, -32.99, 8.44),
         ('SC', 48.625, -81.454, 0.6516), ('SL', 63.121, -97.084, 0.8558), ('SO', 58.789, -81.637, 0.9783),
         ('SS', 19.529, 135.017, 1.478), ('SU', 55.442, -72.955, 1.103), ('SW', -19.019, -39.64, 1.84),
         ('TI', 19.524, 112.175, 1.514), ('TO', 28.807, 2.263, 9.3), ('WL', 22.134, 132.33, 1.546),
         ('YA', 69.067, -97.718, 0.9983))
SLABLESS = ('AF', 'AN', 'NA', 'EU', 'SA')
SLAB_PLATES = ('PA', 'NZ', 'CO', 'AU', 'PS')


@njit
def _crossings(p, v, q):
    """Parity of the crossings of the minor arc from each point of p (m x 3) to q with the closed ring v (k x 3)."""
    out = np.zeros(p.shape[0], np.bool_)
    k = v.shape[0]
    for i in range(p.shape[0]):
        px, py, pz = p[i, 0], p[i, 1], p[i, 2]
        n0 = py*q[2]-pz*q[1]; n1 = pz*q[0]-px*q[2]; n2 = px*q[1]-py*q[0]
        odd = False
        for e in range(k):
            a = v[e]; b = v[(e+1) % k]
            sa = n0*a[0]+n1*a[1]+n2*a[2]; sb = n0*b[0]+n1*b[1]+n2*b[2]
            if (sa > 0) == (sb > 0):
                continue
            m0 = a[1]*b[2]-a[2]*b[1]; m1 = a[2]*b[0]-a[0]*b[2]; m2 = a[0]*b[1]-a[1]*b[0]
            sp = m0*px+m1*py+m2*pz; sq = m0*q[0]+m1*q[1]+m2*q[2]
            if (sp > 0) == (sq > 0):
                continue
            # the great circles meet at +-c; the arcs cross if one of them lies on both
            c0 = n1*m2-n2*m1; c1 = n2*m0-n0*m2; c2 = n0*m1-n1*m0
            for sg in (1., -1.):
                x0, x1, x2 = sg*c0, sg*c1, sg*c2
                if ((x0*(px+q[0])+x1*(py+q[1])+x2*(pz+q[2])) > 0 and
                        (x0*(a[0]+b[0])+x1*(a[1]+b[1])+x2*(a[2]+b[2])) > 0):
                    odd = not odd
                    break
        out[i] = odd
    return out


def _rings():
    js = json.loads((DATA/'tectonicplates-master'/'GeoJSON'/'PB2002_plates.json').read_text())
    rings = {}
    for f in js['features']:
        gm = f['geometry']
        polys = gm['coordinates'] if gm['type'] == 'MultiPolygon' else [gm['coordinates']]
        for poly in polys:
            r = np.asarray(poly[0], float)
            if np.allclose(r[0], r[-1]):
                r = r[:-1]
            v = xyz_from_lonlat(r[:, 0], r[:, 1])
            keep = np.r_[True, np.linalg.norm(np.diff(v, axis=0), axis=1) > 1e-12]
            rings.setdefault(f['properties']['Code'], []).append(np.ascontiguousarray(v[keep]))
    return rings


def _inside(points, v):
    """Points inside the ring v (each part of Bird's plates is smaller than a hemisphere): crossings of the arc to a
    point outside it, the antipode of the ring's mean vertex nudged off any symmetry."""
    c = v.mean(0); c /= np.linalg.norm(c)
    q = -c+np.array([1e-3, -2e-3, 1.5e-3]); q /= np.linalg.norm(q)
    rad = np.arccos(np.clip(v@c, -1, 1)).max()
    near = np.flatnonzero(points@c > np.cos(min(rad+.01, np.pi)))
    out = np.zeros(len(points), bool)
    if near.size:
        out[near] = _crossings(np.ascontiguousarray(points[near]), v, q)
    return out


def plate_of(points, rings, names):
    """Index into names of the plate containing each point (-1: none; the last plate wins an overlap)."""
    best = np.full(len(points), -1)
    for k, nm in enumerate(names):
        hit = np.zeros(len(points), bool)
        for v in rings[nm]:
            hit ^= _inside(points, v)
        best[hit] = k
    return best


def earth_plates(n=217, refresh=False):
    """Bird's plates on CubeSphere(n): (graph, plate index per node, names, continental share per node)."""
    g = CubeSphere(n)
    names = [p[0] for p in POLES]
    cache = OUTPUTS/'planet-reference'/f'earth-plates-n{n}.npz'
    if cache.exists() and not refresh:
        z = np.load(cache)
        pid = z['plate']
    else:
        pid = plate_of(g.xyz, _rings(), names)
        nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
        while (pid < 0).any():                    # gaps along the polygons' seams: the neighbouring plate
            m = pid < 0
            cand = pid[nb[m]].max(1)
            pid[np.flatnonzero(m)] = cand
        np.savez_compressed(cache, plate=pid)
    ref = OUTPUTS/'planet-reference'/'earth-n217.npz'
    if n == 217 and ref.exists():
        cont = (np.load(ref)['elevation_m'] > -1000).astype(float)
    else:
        cont = np.zeros(g.size)
    return g, pid, names, cont


def earth_omega(g, pid, names):
    """Bird's poles as angular velocity (rad/s) in the no-net-rotation frame: {index: vector}."""
    w = {}
    for k, (nm, la, lo, rate) in enumerate(POLES):
        w[k] = xyz_from_lonlat(np.array([lo]), np.array([la]))[0]*np.radians(rate)/1e6/YR
    return remove_net_rotation(g, pid, w)


def remove_net_rotation(g, pid, w):
    a, x = g.area, g.xyz
    I = np.zeros((3, 3)); b = np.zeros(3)
    for k in np.unique(pid):
        m = pid == k
        Ik = a[m].sum()*np.eye(3)-np.einsum('i,ij,ik->jk', a[m], x[m], x[m])
        I += Ik; b += Ik@w.get(int(k), np.zeros(3))
    w0 = np.linalg.solve(I, b)
    return {k: v-w0 for k, v in w.items()}


def earth_steps(names, rings):
    """Bird's steps: midpoint, length (m), class, the two plates (index into names) and the normal out of the first,
    the downgoing plate of a SUB step, Bird's divergence rate (mm/yr, + opening)."""
    s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
    a = xyz_from_lonlat(s['STARTLONG'], s['STARTLAT']); b = xyz_from_lonlat(s['FINALLONG'], s['FINALLAT'])
    m = a+b; m /= np.linalg.norm(m, axis=1, keepdims=True)
    t = b-a; t -= (t*m).sum(1, keepdims=True)*m; t /= np.linalg.norm(t, axis=1, keepdims=True)
    nrm = np.cross(m, t)
    ix = {nm: k for k, nm in enumerate(names)}
    A = np.zeros(len(m), int); B = np.zeros(len(m), int); down = np.full(len(m), -1)
    for i, pb in enumerate(s['PLATEBOUND']):
        for sep in ('/', '\\', '-'):
            if sep in pb:
                x, y = pb.split(sep)
                break
        A[i], B[i] = ix[x], ix[y]
        if sep == '/':
            down[i] = B[i]
        elif sep == '\\':
            down[i] = A[i]
    # which side is A: winding of A's rings at a point 25 km along +nrm
    q = m+25e3/R*nrm; q /= np.linalg.norm(q, axis=1, keepdims=True)
    side = np.zeros(len(m), bool)
    for k in np.unique(A):
        sel = np.flatnonzero(A == k)
        h = np.zeros(sel.size, bool)
        for v in rings[names[k]]:
            h ^= _inside(q[sel], v)
        side[sel] = h
    nA = np.where(side[:, None], -nrm, nrm)                   # out of A
    return {'x': m, 'L': s['STEPLENGTH']*1e3, 'cls': s['STEPCLASS'], 'A': A, 'B': B, 'nA': nA, 'down': down,
            'div': s['VELOCITYDI']}


def terms(g, pid, cont, st, K):
    """The balance's per-plate vectors and matrices (per unit coefficient)."""
    S = np.zeros((K, 3)); P = np.zeros((K, 3)); U = np.zeros((K, 3))
    D = np.zeros((K, 3, 3)); Kc = np.zeros((K, 3, 3)); Z = np.zeros((K, 3, 3))
    X = []                                                     # (i, j, 3x3)
    for k in range(K):
        m = pid == k
        xx = np.einsum('i,ij,ik->jk', g.area[m], g.xyz[m], g.xyz[m])
        D[k] = R**2*(g.area[m].sum()*np.eye(3)-xx)                      # m^4 (graph areas in m^2)
        xc = np.einsum('i,ij,ik->jk', g.area[m]*cont[m], g.xyz[m], g.xyz[m])
        Kc[k] = R**2*((g.area[m]*cont[m]).sum()*np.eye(3)-xc)
    for i in range(len(st['L'])):
        a, b, L, x, nA, c = st['A'][i], st['B'][i], st['L'][i], st['x'][i], st['nA'][i], st['cls'][i]
        qA = np.cross(x, nA)                                  # torque arm of a unit force out of A, per R
        if c == 'SUB':
            d = st['down'][i]; o = b if d == a else a
            qd = qA if d == a else -qA                         # out of the downgoing plate
            S[d] += L*R*qd
            U[o] += L*R*(-qd)                                  # the overriding plate's outward normal is -n_down
            Z[d] += L*R**2*np.outer(qd, qd)
        elif c == 'OSR':
            P[a] += L*R*(-qA); P[b] += L*R*qA
        elif c in ('CCB', 'OCB'):
            X.append((a, b, L*R**2*np.outer(qA, qA)))
    return S, P, U, D, Kc, Z, X


def assemble(T, K):
    """The balance's matrices once: block-diagonal basal (D), keel (K), slab resistance (Z), the coupling (X), and the
    three force vectors, all (3K x 3K) or (3K)."""
    S, P, U, D, Kc, Z, X = T
    def bd(B):
        M = np.zeros((3*K, 3*K))
        for k in range(K):
            M[3*k:3*k+3, 3*k:3*k+3] = B[k]
        return M
    Xm = np.zeros((3*K, 3*K))
    for a, b, M in X:
        Xm[3*a:3*a+3, 3*a:3*a+3] += M; Xm[3*b:3*b+3, 3*b:3*b+3] += M
        Xm[3*a:3*a+3, 3*b:3*b+3] -= M; Xm[3*b:3*b+3, 3*a:3*a+3] -= M
    return {'D': bd(D), 'K': bd(Kc), 'Z': bd(Z), 'X': Xm, 'S': S.reshape(-1), 'P': P.reshape(-1), 'U': U.reshape(-1)}


def solve(Am, prm, K):
    """Angular velocities (rad/s) for coefficients prm = (C_b, beta, F_sp, C_s, F_su, C_c)."""
    Cb, beta, Fsp, Cs, Fsu, Cc = prm
    A = Cb*(Am['D']+beta*Am['K'])+Cs*Am['Z']+Cc*Am['X']
    w = np.linalg.solve(A, Fsp*Am['S']+F_RP*Am['P']+Fsu*Am['U'])
    return {k: w[3*k:3*k+3] for k in range(K)}


class Metric:
    """Per-plate area tensors sum a (I - x x^T): the integral of |w x x|^2 over a plate is w^T M w; the net rotation
    uses the same tensors."""
    def __init__(self, g, pid, K):
        self.M = np.zeros((K, 3, 3))
        for k in range(K):
            m = pid == k
            self.M[k] = g.area[m].sum()*np.eye(3)-np.einsum('i,ij,ik->jk', g.area[m], g.xyz[m], g.xyz[m])
        self.I = self.M.sum(0)
        self.K = K

    def nnr(self, w):
        w0 = np.linalg.solve(self.I, sum(self.M[k]@w[k] for k in range(self.K)))
        return {k: w[k]-w0 for k in range(self.K)}

    def misfit(self, w, wo):
        num = sum((w[k]-wo[k])@self.M[k]@(w[k]-wo[k]) for k in range(self.K))
        return num/sum(wo[k]@self.M[k]@wo[k] for k in range(self.K))


def velocity_misfit(g, pid, w, w_obs):
    return Metric(g, pid, len(w_obs)).misfit(w, w_obs)


def rms_speed(g, pid, w, k):
    m = pid == k
    v = np.cross(w[k], g.xyz[m])*R*YR*100
    return float(np.sqrt(np.average((v*v).sum(1), weights=g.area[m])))


def step_rates(st, w):
    """Divergence rate (mm/yr, + opening) of each step from the poles w: B relative to A along the normal out of A."""
    wa = np.stack([w[int(a)] for a in st['A']]); wb = np.stack([w[int(b)] for b in st['B']])
    v = np.cross(wb-wa, st['x'])*R*YR*1e3
    return (v*st['nA']).sum(1)           # nA points out of A toward B: B moving along it opens the boundary


def wp50(v, L):
    o = np.argsort(v); c = np.cumsum(L[o])
    return float(v[o][np.searchsorted(c, c[-1]/2)])


def nelder_mead(f, x0, step, iters=600):
    n = len(x0)
    pts = [np.asarray(x0, float)]+[np.asarray(x0, float)+np.eye(n)[i]*step[i] for i in range(n)]
    val = [f(p) for p in pts]
    for _ in range(iters):
        o = np.argsort(val); pts = [pts[i] for i in o]; val = [val[i] for i in o]
        c = np.mean(pts[:-1], 0)
        xr = c+(c-pts[-1]); fr = f(xr)
        if fr < val[0]:
            xe = c+2*(c-pts[-1]); fe = f(xe)
            pts[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            pts[-1], val[-1] = xr, fr
        else:
            xc = c+.5*(pts[-1]-c); fc = f(xc)
            if fc < val[-1]:
                pts[-1], val[-1] = xc, fc
            else:
                pts = [pts[0]]+[pts[0]+.5*(p-pts[0]) for p in pts[1:]]; val = [val[0]]+[f(p) for p in pts[1:]]
    i = int(np.argmin(val))
    return pts[i], val[i]


# the fit's box (literature ranges; the fit may not leave them): log10 C_b (Pa s/m), log10 beta, log10 F_sp (N/m),
# log10 terminal slab speed F_sp/C_s (cm/yr), suction as a fraction of F_sp, log10 C_c / C_s
BOX = ((13., 15.5), (-1., 1.5), (12.7, 13.7), (np.log10(2.), np.log10(30.)), (0., 1.), (-3., 1.))


def unpack(z, suction, coupling):
    """Fit coordinates (clipped to BOX) -> (C_b, beta, F_sp, C_s, F_su, C_c)."""
    z = np.array([min(max(v, lo), hi) for v, (lo, hi) in zip(z, BOX)])
    Cb, beta, Fsp = 10**z[0], 10**z[1], 10**z[2]
    Cs = Fsp/(10**z[3]*1e-2/YR)
    Fsu = z[4]*Fsp if suction else 0.
    Cc = 10**z[5]*Cs if coupling else 0.
    return Cb, beta, Fsp, Cs, Fsu, Cc


def box_penalty(z):
    return sum(max(lo-v, 0.)+max(v-hi, 0.) for v, (lo, hi) in zip(z, BOX))


def length_law_omega(g, pid, cont, st, K, p=PARAMS):
    """The current law ('lengths'), fully relaxed, on Earth's plates: speed from Bird's lengths, pole from the
    normalised pull (slab nodes toward trenches, ridge_torque x ridge nodes away)."""
    w = {}
    Ls = np.zeros(K); Lr = np.zeros(K); Lc = np.zeros(K); pull = np.zeros((K, 3))
    for i in range(len(st['L'])):
        a, b, L, x, c = st['A'][i], st['B'][i], st['L'][i]/1e3, st['x'][i], st['cls'][i]
        if c == 'SUB':
            d = st['down'][i]; Ls[d] += L; pull[d] += L*np.cross(_centroid(g, pid, d), x)
        elif c in ('OSR', 'CRB'):
            Lr[a] += L; Lr[b] += L
            for k in (a, b):
                pull[k] += p['ridge_torque']*L*np.cross(x, _centroid(g, pid, k))
        elif c == 'CCB':
            Lc[a] += L; Lc[b] += L
    for k in range(K):
        m = pid == k
        cf = np.average(cont[m], weights=g.area[m]) if m.any() else 0.
        t = (p['v_min_cm_yr']+(p['v_slab_cm_yr']-p['v_min_cm_yr'])*min(1., Ls[k]/p['slab_ref_km'])
             + p['v_ridge_cm_yr']*min(1., Lr[k]/p['ridge_ref_km']))
        t *= (1-p['cont_drag']*cf)/(1+Lc[k]/p['collision_ref_km'])
        ax = pull[k]/np.linalg.norm(pull[k]) if np.linalg.norm(pull[k]) > 0 else np.array([0., 0., 1.])
        w[k] = ax*t*1e-2/R/YR
    return w


_CENT = {}


def _centroid(g, pid, k):
    if k not in _CENT:
        m = pid == k
        c = (g.area[m, None]*g.xyz[m]).sum(0); _CENT[k] = c/np.linalg.norm(c)
    return _CENT[k]


def score(g, pid, names, st, w, w_obs):
    ix = {nm: k for k, nm in enumerate(names)}
    vr = 1-velocity_misfit(g, pid, w, w_obs)
    sp = {nm: (rms_speed(g, pid, w, ix[nm]), rms_speed(g, pid, w_obs, ix[nm])) for nm in names}
    r = step_rates(st, w)
    cls, L = st['cls'], st['L']
    rates = {}
    for c, sg in (('OSR', 1), ('SUB', -1), ('OCB', -1), ('CCB', -1)):
        k = cls == c
        rates[c] = (wp50(sg*r[k], L[k]), wp50(sg*st['div'][k], L[k]))
    ok = {'T-a': vr >= .6,
          'T-b': all(sp[nm][0] <= 3.5 for nm in SLABLESS),
          'T-c': all(1/1.35 <= sp[nm][0]/sp[nm][1] <= 1.35 for nm in SLAB_PLATES),
          'T-d': (abs(rates['OSR'][0]/rates['OSR'][1]-1) <= .25 and abs(rates['SUB'][0]/rates['SUB'][1]-1) <= .25
                  and rates['OCB'][0] <= 20 and .5 <= rates['CCB'][0]/rates['CCB'][1] <= 2)}
    return {'vr': vr, 'speeds': sp, 'rates': rates, 'pass': ok}


def fit(n=217, suction=True, coupling=True, verbose=True):
    g, pid, names, cont = earth_plates(n)
    K = len(names)
    w_obs = earth_omega(g, pid, names)
    st = earth_steps(names, _rings())
    T = terms(g, pid, cont, st, K)
    Am = assemble(T, K)
    met = Metric(g, pid, K)

    def f(z):
        try:
            w = met.nnr(solve(Am, unpack(z, suction, coupling), K))
        except np.linalg.LinAlgError:
            return 1e9
        return met.misfit(w, w_obs)+box_penalty(z)
    z0 = np.array([15., .5, 13.5, 1., .2, -1.])
    best = None
    for start in (z0, z0+[-.5, -.5, -.3, -.3, .3, 1.], z0+[.3, .5, 0, .3, 0, -1.]):
        z, v = nelder_mead(f, start, [.3, .3, .3, .2, .2, .5])
        z, v = nelder_mead(f, z, [.1]*6)
        if best is None or v < best[1]:
            best = (z, v)
    z = best[0]
    prm = unpack(z, suction, coupling)
    w = met.nnr(solve(Am, prm, K))
    return {'Am': Am, 'met': met, 'g': g, 'pid': pid, 'names': names, 'cont': cont, 'st': st, 'T': T, 'w_obs': w_obs, 'w': w, 'prm': prm,
            'z': z, 'K': K}


def params_coefficients(p=PARAMS):
    """planet_tectonics' torque-law PARAMS as (C_b, beta, F_sp, C_s, F_su, C_c)."""
    Cs = p['slab_pull_n_m']/(p['slab_terminal_cm_yr']*1e-2/YR)
    return (p['basal_drag_pa_s_m'], p['keel_drag_factor']-1, p['slab_pull_n_m'], Cs,
            p['slab_suction_frac']*p['slab_pull_n_m'], p['convergent_coupling']*Cs)


def evaluate_params(n=217, p=PARAMS):
    """The torque law with planet_tectonics' PARAMS on Earth's plates: the fit's result dict, coefficients fixed."""
    assert p['ridge_push_n_m'] == F_RP
    g, pid, names, cont = earth_plates(n)
    K = len(names)
    st = earth_steps(names, _rings())
    T = terms(g, pid, cont, st, K)
    Am, met = assemble(T, K), Metric(g, pid, K)
    prm = params_coefficients(p)
    return {'g': g, 'pid': pid, 'names': names, 'cont': cont, 'st': st, 'T': T, 'w_obs': earth_omega(g, pid, names),
            'w': met.nnr(solve(Am, prm, K)), 'prm': prm, 'K': K, 'Am': Am, 'met': met}


def report(r, label):
    g, pid, names, st, w, w_obs = r['g'], r['pid'], r['names'], r['st'], r['w'], r['w_obs']
    sc = score(g, pid, names, st, w, w_obs)
    print(f'\n{label}')
    if 'prm' in r:
        Cb, beta, Fsp, Cs, Fsu, Cc = r['prm']
        print(f'  C_b {Cb:.3g} Pa s/m (eta/h; 1e20 Pa s / 150 km = 6.7e14), continental keel x(1+beta) = x{1+beta:.2f};'
              f' F_sp {Fsp:.3g} N/m, F_rp {F_RP:.3g} (fixed), F_su {Fsu:.3g}; C_s {Cs:.3g} N s/m^2'
              f' (terminal slab speed F_sp/C_s {Fsp/Cs*YR*100 if Cs else float("nan"):.1f} cm/yr); C_c {Cc:.3g} N s/m^2')
    print(f"  T-a variance reduction {sc['vr']:.3f}  {'pass' if sc['pass']['T-a'] else 'FAIL'}")
    print('  T-b slab-less: ' + ', '.join(f'{nm} {sc["speeds"][nm][0]:.1f} (Earth {sc["speeds"][nm][1]:.1f})' for nm in SLABLESS)
          + f"  {'pass' if sc['pass']['T-b'] else 'FAIL'}")
    print('  T-c slab plates: ' + ', '.join(f'{nm} {sc["speeds"][nm][0]:.1f} ({sc["speeds"][nm][1]:.1f})' for nm in SLAB_PLATES)
          + f"  {'pass' if sc['pass']['T-c'] else 'FAIL'}")
    print('  T-d rates p50 law (Bird), mm/yr: ' + ', '.join(f'{c} {a:.0f} ({b:.0f})' for c, (a, b) in sc['rates'].items())
          + f"  {'pass' if sc['pass']['T-d'] else 'FAIL'}")
    big = sorted(range(r['K']), key=lambda k: -(g.area[pid == k].sum()))[:16]
    print('  plates by area, law (Earth) cm/yr: ' + ', '.join(f'{names[k]} {sc["speeds"][names[k]][0]:.1f} ({sc["speeds"][names[k]][1]:.1f})' for k in big))
    return sc


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--n', type=int, default=217)
    ap.add_argument('--json')
    a = ap.parse_args(argv)
    out = {}
    base = None
    for suction, coupling in ((True, True), (False, True), (True, False), (False, False)):
        r = fit(a.n, suction, coupling)
        lab = f'torque law, suction {suction}, convergent coupling {coupling}'
        sc = report(r, lab)
        out[lab] = {'prm': r['prm'], 'score': sc}
        base = r
    r = evaluate_params(a.n)
    out['PARAMS torque'] = {'prm': r['prm'], 'score': report(r, "planet_tectonics PARAMS, force_law 'torque' (rounded fit)")}
    g, pid, names, cont, st, K = base['g'], base['pid'], base['names'], base['cont'], base['st'], base['K']
    wl = remove_net_rotation(g, pid, length_law_omega(g, pid, cont, st, K))
    sc = report({k: v for k, v in {**base, 'w': wl}.items() if k != 'prm'},
                "current length law ('lengths', PARAMS), fully relaxed, Bird's lengths")
    out['lengths'] = {'score': sc}
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)))
    return out


if __name__ == '__main__':
    main()
