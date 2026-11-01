"""D (journal 66, forty-fifth pass): an instantaneous mantle flow driven by the planet's own slab record.

Three parts:

1. **Flow kernels.** Incompressible Stokes flow in a spherical shell whose viscosity varies only with radius, solved
   per spherical-harmonic degree l by propagator matrices (Hager & O'Connell 1981). The poloidal state vector is
   X = (u_r, u_h, r sigma_rr / eta0, r sigma_rtheta / eta0) with the velocity u_r Y r^ + u_h grad_1 Y; it obeys
   dX / d ln r = A(l, eta / eta0) X, and inside a layer of constant viscosity the propagator is the matrix exponential,
   evaluated through the eigen-decomposition of A (eigenvalues l+1, l-1, -l, -l-2: the Stokes solutions). A density
   sheet sigma (kg/m^2) at radius r_i is a jump of r_i g sigma / eta0 in the third component. Free slip at the surface
   and at the core-mantle boundary; uniform g; no self-gravitation (a ~10% effect at the lowest degrees). For a unit
   sheet at each load depth the kernels give the surface normal stress (as air-loaded dynamic topography,
   h = -sigma_rr / (rho_m g): the radial tension under a surface held flat pulls a free surface down; a sheet at the
   surface is compensated, h = -sigma / rho_m) and the free-slip surface velocity coefficient. Viscosity: a radial profile file
   (`visc.sh08` in data/mantle/: non-dimensional radius of each layer's base, Pa s).

2. **Slab record** (mantle frame). Each step every oceanic trench node of the downgoing plate puts a parcel into the
   mantle: the area consumed (closing rate x dt x BOUNDARY_NODE_H h) times the lithosphere's thermal mass excess, which
   by isostasy is the floor's own subsidence (rho_m - rho_w)(GDH1(age) - GDH1(0)). In the upper mantle a parcel moves
   along the mean deep slab dip (58 deg; Lallemand et al. 2005) at its closing rate, toward the overriding side; below
   660 km it sinks vertically at 12 mm/yr (van der Meer et al. 2010, 2018) and leaves the record at 2,850 km.
   Lower-mantle parcels merge by (nearest node, load layer), conserving mass, first moment of depth and mean birth time.

3. **Flow.** The parcels below the cut (300 km: the slab's top is the plate's own; below cratonic keels) are summed per
   50 km load layer into spherical-harmonic surface densities (adjoint transform), and the kernels give the air-loaded
   dynamic topography and the free-slip surface velocity at the nodes, degrees 1..lmax (31: the Earth models'
   truncation).

Nothing here is fitted: the viscosity and the slab kinematics are Earth measurements, the mass is the planet's own
age-depth curve.
"""
import numpy as np

from planet_grid import fields

G = 9.81
GRAV_CONST = 6.674e-11            # m^3 kg^-1 s^-2 (the geoid kernels; pass 50)
CMB_RADIUS_M = 3480e3
RHO_CORE_JUMP = 4340.             # kg/m^3: PREM's density step at the CMB (9,903 - 5,566; Dziewonski & Anderson 1981),
                                  # used only by the self-gravitating geoid bound
ETA0 = 1e21                       # scaling viscosity of the non-dimensional system (any value; results do not depend on it)
Z_660_M = 660e3                   # upper / lower mantle boundary
Z_EXIT_M = 2850e3                 # parcels leave the record here (the slab graveyard above the CMB)
LAYER_M = 50e3                    # load layers below the cut


# ------------------------------------------------------------------ spherical-harmonic helpers (ducc0, lmax-major)
def alm_index(l, m, lmax):
    return m*(2*lmax+1-m)//2+l


def alm_l(lmax):
    """Degree l of every entry of a ducc0 alm array (m-major, m <= l)."""
    return np.concatenate([np.arange(m, lmax+1) for m in range(lmax+1)])


def alm_m(lmax):
    return np.concatenate([np.full(lmax+1-m, m) for m in range(lmax+1)])


def spectrum(alm, lmax):
    """Mean-square power per degree, P_l = (|a_l0|^2 + 2 sum_m>0 |a_lm|^2) / 4 pi (orthonormal harmonics), so
    sum_l P_l is the field's mean square over the sphere."""
    l, m = alm_l(lmax), alm_m(lmax)
    w = np.where(m == 0, 1., 2.)
    return np.bincount(l, w*np.abs(alm)**2, lmax+1)/(4*np.pi)


def gradient_power(alm, lmax):
    """Mean square of the surface gradient (unit sphere) of the field: sum_l l(l+1) P_l."""
    P = spectrum(alm, lmax)
    l = np.arange(lmax+1)
    return float((l*(l+1.)*P).sum())


def synthesis(alm, lmax, xyz, nthreads=1):
    import ducc0
    return ducc0.sht.experimental.synthesis_general(alm=np.asarray(alm, np.complex128)[None], spin=0, lmax=lmax,
                                                    loc=fields.node_loc(xyz), epsilon=1e-10, nthreads=nthreads)[0]


def gradient(alm, lmax, xyz, nthreads=1):
    """Surface gradient (on the unit sphere) of the scalar with coefficients alm, as xyz vectors at the points: a spin-1
    synthesis of the gradient (E-mode) coefficients sqrt(l(l+1)) a_lm."""
    import ducc0
    loc = fields.node_loc(xyz)
    a = np.asarray(alm, np.complex128)
    e = a*np.sqrt(alm_l(lmax)*(alm_l(lmax)+1.))
    gt, gp = ducc0.sht.experimental.synthesis_general(alm=np.stack([e, 0*e]), spin=1, lmax=lmax, loc=loc,
                                                      epsilon=1e-10, nthreads=nthreads)
    th, ph = loc[:, 0], loc[:, 1]
    eth = np.stack([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), -np.sin(th)], 1)
    eph = np.stack([-np.sin(ph), np.cos(ph), 0*ph], 1)
    return gt[:, None]*eth+gp[:, None]*eph


def point_alm(xyz, weight, lmax, nthreads=1):
    """sum_i weight_i conj(Y_lm(x_i)): the coefficients of a field made of delta functions (per steradian)."""
    import ducc0
    if len(weight) == 0:
        return np.zeros((lmax+1)*(lmax+2)//2, np.complex128)
    return ducc0.sht.experimental.adjoint_synthesis_general(map=np.asarray(weight, np.float64)[None], spin=0, lmax=lmax,
                                                            loc=fields.node_loc(xyz), epsilon=1e-10,
                                                            nthreads=nthreads)[0]


# ------------------------------------------------------------------ viscosity and propagator kernels
def read_viscosity(path):
    """A radial viscosity profile in the HC format: rows (non-dimensional radius of the layer's base, viscosity Pa s),
    bottom up. Returns (base radii as fractions of the planet radius, viscosities)."""
    v = np.loadtxt(path)
    return v[:, 0], v[:, 1]


def system_matrix(l, e):
    """A(l, eta/eta0) of dX / d ln r = A X (Hager & O'Connell 1981; X = (u_r, u_h, r s_rr/eta0, r s_rt/eta0))."""
    L = l*(l+1.)
    return np.array([[-2., L, 0., 0.], [-1., 1., 0., 1./e], [12.*e, -6.*L*e, 1., L], [-6.*e, 2.*(2.*L-1.)*e, -1., -2.]])


def propagator_eig(l, e, x):
    """exp(A x) through the eigen-decomposition of A (real eigenvalues l+1, l-1, -l, -l-2)."""
    w, W = np.linalg.eig(system_matrix(l, e))
    w, W = w.real, W.real
    return (W*np.exp(w*x))@np.linalg.inv(W)


def propagator_rk4(l, e, x, n=4000):
    """exp(A x) by fourth-order Runge-Kutta (the unit test's reference)."""
    A = system_matrix(l, e)
    P = np.eye(4); h = x/n
    for _ in range(n):
        k1 = A@P; k2 = A@(P+h/2*k1); k3 = A@(P+h/2*k2); k4 = A@(P+h*k3)
        P = P+h/6*(k1+2*k2+2*k3+k4)
    return P


class FlowKernels:
    """Surface response to a unit density sheet (1 kg/m^2) at each load radius, per degree l = 1..lmax:
    KH[l, j] = -sigma_rr(surface) / g (kg/m^2; the air-loaded topography is KH sigma / rho_m: negative over a dense
    sheet, -> -1 for a sheet at the surface, -> 0 at the CMB) and
    KV[l, j] = u_h(surface) (m/s per kg/m^2; horizontal velocity u_h grad_1 Y), and (pass 50, for the geoid)
    KC[l, j] = +sigma_rr(CMB) / g (kg/m^2): the mass per area of the CMB's deflection (the core is a hydrostatic fluid, so
    its pressure at the CMB radius has no l >= 1 part and the boundary moves by sigma_rr / (drho_c g); a dense sheet
    pushes it down, -> -1 for a sheet at the CMB, -> 0 at the surface)."""

    def __init__(self, radius_m, load_radii_m, base_frac, eta, lmax, cmb_m=CMB_RADIUS_M):
        self.a, self.c, self.lmax = float(radius_m), float(cmb_m), int(lmax)
        tops = np.append(base_frac[1:]*radius_m, radius_m)
        self.tops, self.eta = tops, np.asarray(eta, np.float64)
        self.r = np.asarray(load_radii_m, np.float64)
        J = len(self.r)
        self.KH = np.zeros((lmax+1, J)); self.KV = np.zeros((lmax+1, J)); self.KC = np.zeros((lmax+1, J))
        for l in range(1, lmax+1):
            dec = []
            for et in self.eta:
                w, W = np.linalg.eig(system_matrix(l, et/ETA0))
                dec.append((w.real, W.real, np.linalg.inv(W.real)))
            Pac = self._between(dec, self.a, self.c)
            M = Pac[np.ix_([0, 3], [1, 2])]
            for j, ri in enumerate(self.r):
                s = self._between(dec, self.a, ri)@np.array([0., 0., ri*G/ETA0, 0.])
                Vc, Rc = np.linalg.solve(M, -s[[0, 3]])
                Xa = Pac@np.array([0., Vc, Rc, 0.])+s
                self.KH[l, j] = -Xa[2]*ETA0/self.a/G
                self.KV[l, j] = Xa[1]
                self.KC[l, j] = Rc*ETA0/self.c/G

    def _between(self, dec, r2, r1):
        """Propagator from r1 up to r2 through the viscosity layers."""
        P = np.eye(4); r = r1
        for top, (w, W, Wi) in zip(self.tops, dec):
            if top <= r:
                continue
            nxt = min(top, r2)
            P = ((W*np.exp(w*np.log(nxt/r)))@Wi)@P
            r = nxt
            if r >= r2:
                break
        return P

    def surface(self, loads):
        """loads (J, nalm): surface density coefficients (kg/m^2) on each load sphere. Returns (H, V): the air-loaded
        -sigma_rr / g coefficients (kg/m^2; divide by rho_m for metres) and the horizontal-velocity potential (m/s)."""
        li = alm_l(self.lmax)
        H = (self.KH[li, :].T*loads).sum(0)
        V = (self.KV[li, :].T*loads).sum(0)
        return H, V

    def geoid_kernels(self, self_grav=False, rho_top=None, rho_core=RHO_CORE_JUMP):
        """(KN, KHs): the surface geoid (m per kg/m^2 of load) and the surface deflection mass (kg/m^2 per kg/m^2),
        per degree and load radius. Without self-gravitation (self_grav False) the deflections are KH and KC and
        N_l = 4 pi G_N a / ((2l+1) g) [ (r/a)^(l+2) + KH + KC (c/a)^(l+2) ] (each sheet's exterior potential).
        With self_grav the deflections are measured from the perturbed equipotentials: m_a = KH sigma + rho_top N(a),
        m_c = KC sigma + rho_core N(c), solved per degree; this is the boundary part of self-gravitation in a mantle
        of uniform density (the reduced pressure p + rho_0 phi absorbs the rest), a bound, not PREM's stratification.
        rho_top: the surface density contrast (kg/m^3; required with self_grav)."""
        a, c = self.a, self.c
        L = np.arange(self.lmax+1)[:, None]
        k = 4*np.pi*GRAV_CONST/((2*L+1.)*G)                       # N = k * (radius x surface density) x geometry
        ra = self.r[None, :]/a
        load_a = self.r[None, :]*ra**(L+1)                       # r_j (r_j/a)^(l+1): sheet seen from the surface
        if not self_grav:
            KN = k*(load_a+a*self.KH+c*self.KC*(c/a)**(L+1))
            KN[0] = 0.
            return KN, self.KH.copy()
        if rho_top is None:
            raise ValueError('self_grav needs rho_top')
        load_c = self.r[None, :]*(c/self.r[None, :])**L          # sheet seen from the CMB (inside it)
        KN = np.zeros_like(self.KH); KHs = np.zeros_like(self.KH)
        for l in range(1, self.lmax+1):
            kl = k[l, 0]; q = (c/a)**(l+1); p = (c/a)**l
            M = np.array([[1-rho_top*kl*a, -rho_top*kl*c*q], [-rho_core*kl*a*p, 1-rho_core*kl*c]])
            rhs = np.stack([self.KH[l]+rho_top*kl*load_a[l], self.KC[l]+rho_core*kl*load_c[l]])
            ma, mc = np.linalg.solve(M, rhs)
            KN[l] = kl*(a*ma+load_a[l]+c*q*mc)
            KHs[l] = ma
        return KN, KHs

    def geoid(self, loads, self_grav=False, rho_top=None, rho_core=RHO_CORE_JUMP):
        """Geoid coefficients (m) of the loads (J, nalm; kg/m^2), degrees 1..lmax (degree 1 is reported but has no
        meaning in a centre-of-mass frame)."""
        KN, _ = self.geoid_kernels(self_grav, rho_top, rho_core)
        li = alm_l(self.lmax)
        return (KN[li, :].T*loads).sum(0)


# ------------------------------------------------------------------ thermal expansivity with depth
ALPHA_GDH1 = 3.1e-5               # K^-1: the expansivity of Stein & Stein's (1992) GDH1 fit, which sets the parcels' mass
ALPHA_UM = ((0., 3.1e-5), (400., 2.378e-5), (660., 2.079e-5))   # km, K^-1: surface (GDH1); Steinberger & Calderwood
                                  # (2006, GJI 167, 1461) phase-boundary averages at 400 and 660 km
ALPHA_LM_TOP, ALPHA_LM_BASE = 2.5e-5, 1.0e-5   # just below 670 km and at the CMB (Steinberger & Calderwood 2006, sect.
                                  # 2.3.5, in agreement with Schubert et al. 2001)
ALPHA_DELTA_T, ALPHA_B = 5.5, 1.4  # their alpha = alpha0 exp[-(delta_T / b)(1 - (rho0 / rho)^b)] (Anderson et al. 1992)


def prem_lower_mantle_density(z_km):
    """PREM density (kg/m^3) for 670-2891 km (Dziewonski & Anderson 1981: 7.9565 - 6.4761 x + 5.5283 x^2 - 3.0807 x^3)."""
    x = (6371.-np.asarray(z_km, np.float64))/6371.
    return 1e3*(7.9565-6.4761*x+5.5283*x**2-3.0807*x**3)


def alpha_ratio(z_km):
    """alpha(z) / ALPHA_GDH1: how much of a slab parcel's surface thermal mass excess remains as density anomaly at
    depth z (the parcel's mass and temperature deficit are conserved; its density anomaly is rho alpha dT per volume,
    and its volume scales as 1 / rho, so the mass anomaly scales as alpha alone). Upper mantle: linear between the
    ALPHA_UM points; lower mantle: Steinberger & Calderwood's form along PREM's density, its two constants (rho0 and
    the prefactor) set by their stated values at the top and the base."""
    z = np.asarray(z_km, np.float64)
    um = np.interp(z, [p[0] for p in ALPHA_UM], [p[1] for p in ALPHA_UM])
    rt, rb = prem_lower_mantle_density(670.), prem_lower_mantle_density(2891.)
    c = ALPHA_DELTA_T/ALPHA_B
    k = np.log(ALPHA_LM_TOP/ALPHA_LM_BASE)/c/(rt**-ALPHA_B-rb**-ALPHA_B)      # rho0^b
    xa = lambda r: k*r**-ALPHA_B
    A = ALPHA_LM_TOP/np.exp(-c*(1-xa(rt)))
    lm = A*np.exp(-c*(1-xa(prem_lower_mantle_density(np.clip(z, 670., 2891.)))))
    return np.where(z < 660., um, np.where(z < 670., np.interp(z, [660., 670.], [ALPHA_UM[-1][1], ALPHA_LM_TOP]), lm))/ALPHA_GDH1


# The slab's thermal mass at depth per unit of its surface (GDH1) mass, from a thermodynamic model fitted to
# experiments (Stixrude & Lithgow-Bertelloni 2011, SLB2011, evaluated with BurnMan 2.1.0; fifty-first pass, journal
# 72): a GDH1 plate (95 km, 1,450 C) of age 63 Myr laid on the ambient isentrope (potential temperature 1,600 K) of a
# pyrolite with fixed modes per region; every element follows its own isentrope (latent-heat jumps 410 +2.5%, 660
# -1.0%); the mass is sum (1 - rho_ambient / rho_element) dm over the elements, normalised by the same sum for the
# plate column at the trench, so the surface mass stays the measured subsidence mass. The mean of the no-conduction
# and fully diffused limits (they differ by <= 2%); potential temperature 1,550-1,650 K and plate ages 42-85 Myr move
# it by <= 2%. Unlike alpha_ratio it keeps the parcel's entropy deficit, not its temperature deficit: along the
# adiabat a parcel's temperature deficit grows with the absolute temperature (x 1.2 at 825 km, x 1.6 at 2,825 km).
# Tabulated at the load-layer mid-depths below a 300 km cut (325, 375, ..., 2,825 km).
THERMO_Z0_KM, THERMO_DZ_KM = 325., 50.
THERMO_F = (0.7106, 0.6756, 0.7352, 0.7115, 0.6594, 0.6403, 0.6222, 0.8388, 0.8245, 0.8110, 0.7980, 0.7857, 0.7739,
            0.7627, 0.7519, 0.7415, 0.7316, 0.7220, 0.7128, 0.7039, 0.6953, 0.6870, 0.6790, 0.6713, 0.6638, 0.6566,
            0.6495, 0.6427, 0.6361, 0.6296, 0.6234, 0.6173, 0.6113, 0.6055, 0.5999, 0.5944, 0.5890, 0.5838, 0.5787,
            0.5737, 0.5688, 0.5640, 0.5593, 0.5546, 0.5501, 0.5457, 0.5413, 0.5371, 0.5328, 0.5287, 0.5246)


def thermo_ratio(z_km):
    """The slab's thermal mass at depth z per unit of its surface mass (THERMO_F, linear between the tabulated
    layer mid-depths, held at the ends)."""
    z = np.asarray(z_km, np.float64)
    return np.interp(z, THERMO_Z0_KM+THERMO_DZ_KM*np.arange(len(THERMO_F)), np.asarray(THERMO_F))


# ------------------------------------------------------------------ the slab record
class SlabRecord:
    """Subducted lithosphere in the mantle frame: parcel position (unit vector), depth (m), mass excess (kg), birth
    time (Myr), and in the upper mantle its along-dip speed (m/yr) and horizontal direction of descent (unit tangent)."""

    def __init__(self, dip_deg=58., lm_sink_mm_yr=12.):
        self.sd, self.cd = np.sin(np.radians(dip_deg)), np.cos(np.radians(dip_deg))
        self.v_lm = lm_sink_mm_yr*1e-3                   # m/yr
        self.pos = np.zeros((0, 3)); self.dirn = np.zeros((0, 3))
        self.z = np.zeros(0); self.mass = np.zeros(0); self.born = np.zeros(0); self.speed = np.zeros(0)
        self.deposited_kg = 0.
        self.exited_kg = 0.

    @property
    def size(self):
        return self.z.size

    def total_mass(self):
        return float(self.mass.sum())

    def deposit(self, xyz, dirn, mass, speed_m_yr, t):
        """New parcels at the surface (z = 0) above the trench."""
        keep = mass > 0
        if not keep.any():
            return
        self.pos = np.concatenate([self.pos, xyz[keep]]); self.dirn = np.concatenate([self.dirn, dirn[keep]])
        self.z = np.concatenate([self.z, np.zeros(keep.sum())]); self.mass = np.concatenate([self.mass, mass[keep]])
        self.born = np.concatenate([self.born, np.full(keep.sum(), float(t))])
        self.speed = np.concatenate([self.speed, speed_m_yr[keep]])
        self.deposited_kg += float(mass[keep].sum())

    def advance(self, dt_myr, radius):
        """Move every parcel over dt: along the dip in the upper mantle (the part of the step left after crossing
        660 km is spent sinking vertically), vertically below; parcels past Z_EXIT_M leave the record."""
        if self.size == 0:
            return
        T = dt_myr*1e6
        um = self.z < Z_660_M
        tu = np.zeros(self.size)
        if um.any():
            vz = self.speed[um]*self.sd
            t_cross = np.where(vz > 0, (Z_660_M-self.z[um])/np.maximum(vz, 1e-30), np.inf)
            tu[um] = np.minimum(T, t_cross)
            dz = vz*tu[um]
            dh = self.speed[um]*self.cd*tu[um]
            th = dh/radius
            p, d = self.pos[um], self.dirn[um]
            c, s = np.cos(th)[:, None], np.sin(th)[:, None]
            self.pos[um] = c*p+s*d
            self.dirn[um] = c*d-s*p
            self.z[um] = np.where(tu[um] >= t_cross, Z_660_M, self.z[um]+dz)
        self.z = self.z+np.where(self.z >= Z_660_M, self.v_lm*(T-tu), 0.)
        out = self.z >= Z_EXIT_M
        if out.any():
            self.exited_kg += float(self.mass[out].sum())
            self._take(~out)

    def _take(self, m):
        self.pos, self.dirn, self.z = self.pos[m], self.dirn[m], self.z[m]
        self.mass, self.born, self.speed = self.mass[m], self.born[m], self.speed[m]

    def merge(self, graph, z_cut):
        """Merge lower-mantle parcels that share the nearest node and the load layer (mass-weighted position, depth
        and birth time)."""
        lm = self.z >= Z_660_M
        if lm.sum() < 2:
            return
        I = np.flatnonzero(lm)
        node = fields.nearest_nodes(graph, self.pos[I])
        lay = np.floor((self.z[I]-z_cut)/LAYER_M).astype(np.int64)
        key = node.astype(np.int64)*4096+lay
        u, inv = np.unique(key, return_inverse=True)
        if u.size == I.size:
            return
        m = self.mass[I]
        M = np.bincount(inv, m, u.size)
        Mw = np.where(M > 0, M, 1.)
        pos = np.stack([np.bincount(inv, m*self.pos[I, k], u.size) for k in range(3)], 1)
        nrm = np.linalg.norm(pos, axis=1)
        first = np.full(u.size, -1, np.int64); first[inv[::-1]] = I[::-1]
        pos = np.where(nrm[:, None] > 0, pos/np.where(nrm > 0, nrm, 1.)[:, None], self.pos[first])
        z = np.bincount(inv, m*self.z[I], u.size)/Mw
        born = np.bincount(inv, m*self.born[I], u.size)/Mw
        z = np.where(M > 0, z, self.z[first]); born = np.where(M > 0, born, self.born[first])
        keep = ~lm
        self.pos = np.concatenate([self.pos[keep], pos]); self.dirn = np.concatenate([self.dirn[keep], self.dirn[first]])
        self.z = np.concatenate([self.z[keep], z]); self.mass = np.concatenate([self.mass[keep], M])
        self.born = np.concatenate([self.born[keep], born]); self.speed = np.concatenate([self.speed[keep], self.speed[first]])

    def layer_loads(self, z_cut, n_layers, radius, lmax, nthreads=1, factor=None):
        """Surface-density coefficients (kg/m^2 on each load sphere) of the parcels below z_cut, per LAYER_M layer;
        factor (per layer) scales each layer's mass (the expansivity with depth)."""
        nalm = (lmax+1)*(lmax+2)//2
        out = np.zeros((n_layers, nalm), np.complex128)
        below = self.z >= z_cut
        if not below.any():
            return out
        lay = np.minimum(((self.z[below]-z_cut)//LAYER_M).astype(np.int64), n_layers-1)
        pos, mass = self.pos[below], self.mass[below]
        for j in np.unique(lay):
            m = lay == j
            r = radius-(z_cut+(j+.5)*LAYER_M)
            out[j] = point_alm(pos[m], mass[m]/r**2*(1. if factor is None else factor[j]), lmax, nthreads)
        return out

    def mass_by_layer(self, z_cut, n_layers):
        below = self.z >= z_cut
        lay = np.minimum(((self.z[below]-z_cut)//LAYER_M).astype(np.int64), n_layers-1)
        return np.bincount(lay, self.mass[below], n_layers)


def load_radii(radius, z_cut):
    """Mid radii of the load layers between z_cut and Z_EXIT_M."""
    n = int(round((Z_EXIT_M-z_cut)/LAYER_M))
    return radius-(z_cut+(np.arange(n)+.5)*LAYER_M)


def thermal_mass_excess(age_myr, rho_mantle, rho_water):
    """Mass excess per area (kg/m^2) of oceanic lithosphere of the given age: the water-replaced mantle of its
    thermal subsidence, (rho_m - rho_w)(GDH1(age) - GDH1(0))."""
    from planet_tectonics.isostasy import gdh1_depth
    return (rho_mantle-rho_water)*(gdh1_depth(age_myr)-gdh1_depth(0.))


class MantleFlow:
    """The record, the kernels and the current flow on a graph. p: rho_mantle, flow_lmax, flow_z_cut_km, slab_dip_deg,
    lm_sink_mm_yr, viscosity_file (path), threads, slab_alpha ('constant': the surface expansivity at all depths;
    'depth': alpha_ratio; 'thermo': thermo_ratio)."""

    def __init__(self, graph, p):
        from planet_grid.paths import DATA
        self.g = graph
        self.lmax = int(p['flow_lmax'])
        self.z_cut = float(p['flow_z_cut_km'])*1e3
        self.rho_m = float(p['rho_mantle'])
        self.nthreads = 1 if p.get('threads') is None else int(p['threads'])
        vf = p.get('viscosity_file') or DATA/'mantle'/'visc.sh08'
        base, eta = read_viscosity(vf)
        self.radii = load_radii(graph.radius, self.z_cut)
        self.kernels = FlowKernels(graph.radius, self.radii, base, eta, self.lmax)
        sa = p.get('slab_alpha', 'thermo')
        if sa == 'thermo':                                # the slab's own thermodynamics (SLB2011, pass 51)
            self.factor = thermo_ratio((graph.radius-self.radii)/1e3)
        elif sa == 'depth':                               # the ambient expansivity's decrease with depth
            self.factor = alpha_ratio((graph.radius-self.radii)/1e3)
        elif sa == 'constant':
            self.factor = None
        else:
            raise ValueError(f"slab_alpha {sa!r}: 'thermo', 'depth' or 'constant'")
        self.record = SlabRecord(p['slab_dip_deg'], p['lm_sink_mm_yr'])
        nalm = (self.lmax+1)*(self.lmax+2)//2
        self.H = np.zeros(nalm, np.complex128)            # air-loaded dynamic topography (m)
        self.V = np.zeros(nalm, np.complex128)            # surface-velocity potential (m/yr on the unit sphere)
        self.h_air = np.zeros(graph.size)
        self.v_xyz = np.zeros((graph.size, 3))            # free-slip surface velocity, m/yr

    def update(self):
        """The flow of the current record."""
        loads = self.record.layer_loads(self.z_cut, len(self.radii), self.g.radius, self.lmax, self.nthreads, self.factor)
        H, V = self.kernels.surface(loads)
        self.H = H/self.rho_m
        self.V = V*3.15576e7
        self.h_air = synthesis(self.H, self.lmax, self.g.xyz, self.nthreads)
        self.v_xyz = gradient(self.V, self.lmax, self.g.xyz, self.nthreads)
