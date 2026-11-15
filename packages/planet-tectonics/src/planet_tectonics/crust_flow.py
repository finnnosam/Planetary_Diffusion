"""The continental crust's lateral flow from a measured rheology (eighty-fourth pass, journal 109).

Today's crust flow ('fitted') diffuses crustal thickness with kappa_normal (30 m^2/yr, cold crust) plus a collapse
term above 45 km, crust_kappa * (exp((T - 45 km) / 6 km) - 1) up to 200 x crust_kappa, whose four numbers were tuned
to hold plateaus near Tibet's and the Altiplano's ~70 km (an M-target, fitted). It holds the planet's belts down
(without it towers, journal 101) but keeps 0.174 of the surface's continental crust above 45 km (Earth 0.086).

This module replaces the collapse term by the flow the crust's own rheology gives:

1. A steady conductive geotherm in each column from its crustal thickness T and lithospheric thickness L (C's
   field): the felsic upper part (1 - lower_crust_share of T) produces flow_heat_upper W/m^3, the lower crust
   flow_heat_lower (Hasterok & Chapman 2011: 1.0 and 0.4 uW/m^3 from exposed sections and granulite terranes; the bulk
   of a 40 km crust is then 0.81, inside Jaupart et al. 2016's 0.65-0.87), none below; the base of the lithosphere is
   at flow_t_lab_c (1,175 C: the isotherm Hoggard et al. 2020 map as L, whose L the planet's frame uses). Thickening
   by pure shear thickens both layers, so a thick column carries more heat sources over a longer path and its base
   runs hot: on a 40 km crust over 178 km of lithosphere the surface heat flux is 47 mW/m^2 and the Moho's 15 (Earth's
   cratons ~40-45 and 10-15, Jaupart & Mareschal 1999); a 60 km crust over 140 km reaches 820 C at its base.
   Steady state is the upper bound (a freshly thickened column heats over ~20-100 Myr).
2. Dislocation creep with laboratory flow laws, felsic part: wet quartzite (Gleason & Tullis 1995, melt-free:
   A 1.1e-4 MPa^-4 s^-1, n 4, Q 223 kJ/mol); lower crust: wet anorthite (Rybacki & Dresen 2000: log A 2.6 MPa^-3 s^-1,
   n 3, Q 356 kJ/mol). Simple shear in the tensor form, gamma_dot = 3^((n+1)/2) A tau^n exp(-Q / R T).
3. Driving stress from the isostatic surface slope: tau(d) = rho_c g S (d - d0), S = (1 - rho_c / rho_m) |grad T|,
   uniform with depth in the crust. Two boundary modes (Rey, Vanderhaeghe & Teyssier 2001):
   'channel' (fixed boundary): the upper crust is held by the plate around it, so the velocity is zero at the
     surface and at the Moho (Poiseuille with depth-dependent viscosity; d0 where the shear integrates to zero);
   'sheet' (free boundary): the whole crust spreads under its own weight on a fixed Moho (d0 = 0, the stress grows
     from the free surface down), the thin-viscous-sheet limit;
   'rey': each mode where Rey et al. put it: fixed while the orogen converges, free once the boundary forces relax.
     An edge is pinned (channel) when either node lies within flow_pinned_km of active convergence (the continental
     arc belt over a converging trench, or a collision boundary closing faster than boundary_min_mm_yr; last step's
     boundaries): 290 km, the distance within which Earth's collisions put half their crustal thickening
     (crust_kernels.COLLISION_KERNEL, GSRM v2.1 on Bird's CCB steps); free (sheet) elsewhere.
   The flux per unit width q = |integral of d gamma_dot(d)| over the crust; the edge diffusivity is q / |grad T|.

The flux is tabulated once in (T, L, log10 S) and interpolated (trilinearly in log q) on every edge of the
continental graph at every substep. kappa_normal is kept (cold-crust creep, Earth estimate). flow_kappa_max_m2_yr is
numerical: at today's collapse ceiling (2e5 m^2/yr) a node-scale step relaxes in ~10 kyr and a 300 km plateau edge in
~0.5 Myr, inside one step; the share of edges at the ceiling is booked (crust_flow_capped_edges). Time stepping is
explicit and conservative; edges whose diffusivity exceeds flow_kappa_split_m2_yr (numerical) take their own substeps
inside each step of the rest (operator splitting), so a few hot plateau edges do not set the step of the whole
continent.
No number here is fitted to the planet."""
import numpy as np

from planet_grid.jit import njit

R_GAS = 8.314
YEAR_S = 3.15576e7
T_AXIS_KM = np.arange(10., 100.5, 1.)
L_AXIS_KM = np.arange(40., 320.5, 5.)
LOGS_AXIS = np.arange(-6., -.49, .25)
NPTS = 200


@njit(cache=True)
def _flux_table(T_ax, L_ax, logS_ax, share, A1, A2, k_c, k_m, t0, tl, law1, law2, rho_c, g, mode, npts):
    """log10 of the flux per unit width (m^2/yr) on the (T km, L km, log10 S) grid."""
    out = np.empty((T_ax.shape[0], L_ax.shape[0], logS_ax.shape[0]))
    d = np.empty(npts)
    c = np.empty(npts)
    nn = np.empty(npts)
    gam = np.empty(npts)
    for it in range(T_ax.shape[0]):
        T = T_ax[it]*1e3
        Tf = T*(1.-share)
        for il in range(L_ax.shape[0]):
            L = max(L_ax[il]*1e3, T+10e3)
            # Moho temperature linear in the Moho heat flux qm: Tm = a + b qm; mantle lid Tm + qm (L - T) / k_m = tl
            q0a = A1*Tf+A2*(T-Tf)
            a = t0+q0a*Tf/k_c-A1*Tf*Tf/(2.*k_c)+(q0a-A1*Tf)*(T-Tf)/k_c-A2*(T-Tf)**2/(2.*k_c)
            b = T/k_c
            qm = (tl-a)/(b+(L-T)/k_m)
            q0 = qm+q0a
            for i in range(npts):
                z = (i+.5)/npts*T
                d[i] = z
                if z < Tf:
                    temp = t0+q0*z/k_c-A1*z*z/(2.*k_c)
                    A, n, Q = law1[0], law1[1], law1[2]
                else:
                    temp = t0+q0*Tf/k_c-A1*Tf*Tf/(2.*k_c)+(q0-A1*Tf)*(z-Tf)/k_c-A2*(z-Tf)**2/(2.*k_c)
                    A, n, Q = law2[0], law2[1], law2[2]
                temp = min(temp, tl)
                c[i] = 3.**((n+1.)/2.)*A*np.exp(-Q/(R_GAS*(temp+273.15)))
                nn[i] = n
            w = T/npts
            for js in range(logS_ax.shape[0]):
                G = rho_c*g*10.**logS_ax[js]/1e6             # MPa per m of depth
                if mode == 1:                                # channel: d0 where the shear integrates to zero
                    lo = 0.
                    hi = T
                    for _ in range(40):
                        mid = .5*(lo+hi)
                        s = 0.
                        for i in range(npts):
                            tau = G*(d[i]-mid)
                            s += c[i]*np.sign(tau)*np.abs(tau)**nn[i]
                        if s > 0.:
                            lo = mid
                        else:
                            hi = mid
                    d0 = .5*(lo+hi)
                else:
                    d0 = 0.
                q = 0.
                for i in range(npts):
                    tau = G*(d[i]-d0)
                    gam[i] = c[i]*np.sign(tau)*np.abs(tau)**nn[i]
                    q += d[i]*gam[i]*w
                q = abs(q)*YEAR_S
                out[it, il, js] = np.log10(max(q, 1e-300))
    return out


def moho_state(T_m, L_m, p):
    """(Moho heat flux, surface heat flux in W/m^2, Moho temperature in C) of the steady geotherm _flux_table uses, for
    crust T_m (m) over lithosphere L_m (m; at least T_m + 10 km)."""
    T = np.asarray(T_m, np.float64)
    L = np.maximum(np.asarray(L_m, np.float64), T+10e3)
    A1, A2, k = p['flow_heat_upper_w_m3'], p['flow_heat_lower_w_m3'], p['flow_k_crust']
    Tf = T*(1.-p['lower_crust_share'])
    q0a = A1*Tf+A2*(T-Tf)
    a = p['flow_t_surface_c']+q0a*Tf/k-A1*Tf**2/(2*k)+(q0a-A1*Tf)*(T-Tf)/k-A2*(T-Tf)**2/(2*k)
    qm = (p['flow_t_lab_c']-a)/(T/k+(L-T)/p['flow_k_mantle'])
    return qm, qm+q0a, a+qm*T/k


_TABLES = {}


def flux_table(p, mode=None):
    """The run's flux table for mode 'channel' or 'sheet' (default: the run's rule; cached by the values it depends
    on)."""
    mode = {'channel': 1, 'sheet': 0}[mode or p['crust_flow_rule']]
    key = (mode, p['lower_crust_share'], p['flow_heat_upper_w_m3'], p['flow_heat_lower_w_m3'], p['flow_k_crust'],
           p['flow_k_mantle'], p['flow_t_surface_c'], p['flow_t_lab_c'], tuple(p['flow_law_upper']),
           tuple(p['flow_law_lower']), p['rho_crust'])
    if key not in _TABLES:
        _TABLES[key] = _flux_table(T_AXIS_KM, L_AXIS_KM, LOGS_AXIS, p['lower_crust_share'], p['flow_heat_upper_w_m3'],
                                   p['flow_heat_lower_w_m3'], p['flow_k_crust'], p['flow_k_mantle'],
                                   p['flow_t_surface_c'], p['flow_t_lab_c'], np.array(p['flow_law_upper'], np.float64),
                                   np.array(p['flow_law_lower'], np.float64), p['rho_crust'], 9.81, mode, NPTS)
    return _TABLES[key]


@njit(cache=True)
def _idx(ax, v):
    if v <= ax[0]:
        return 0, 0.
    if v >= ax[-1]:
        return ax.shape[0]-2, 1.
    h = ax[1]-ax[0]
    i = int((v-ax[0])/h)
    if i > ax.shape[0]-2:
        i = ax.shape[0]-2
    return i, (v-ax[i])/h


@njit(cache=True)
def _interp(tab, T_ax, L_ax, S_ax, T, L, logS):
    """Trilinear interpolation of tab, clamped to the axes."""
    i, fi = _idx(T_ax, T)
    j, fj = _idx(L_ax, L)
    k, fk = _idx(S_ax, logS)
    v = 0.
    for a in range(2):
        wa = fi if a else 1.-fi
        for b in range(2):
            wb = fj if b else 1.-fj
            for cc in range(2):
                wc = fk if cc else 1.-fk
                v += wa*wb*wc*tab[i+a, j+b, k+cc]
    return v


@njit(cache=True)
def _edge_k(ti, tj, d, li, lj, tb, T_ax, L_ax, S_ax, airy, smin):
    """The flow part of one edge's diffusivity, q / |grad T| (m^2/yr), uncapped."""
    grad = abs(tj-ti)/d
    S = airy*grad
    Tm = .5*(ti+tj)/1e3
    Lm = .5*(li+lj)
    if S > smin:
        return 10.**_interp(tb, T_ax, L_ax, S_ax, Tm, Lm, np.log10(S))/grad
    # below the table: q ~ S^n, n >= 3, so kappa falls as S^2 or faster
    return 10.**_interp(tb, T_ax, L_ax, S_ax, Tm, Lm, S_ax[0])/(smin/airy)*(S/smin)**2


@njit(cache=True)
def reverse_index(nbr):
    """rev[i, k]: the slot of i in its neighbour j = nbr[i, k]'s row (-1 where there is none)."""
    rev = np.full(nbr.shape, -1, np.int64)
    for i in range(nbr.shape[0]):
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0:
                continue
            for kk in range(nbr.shape[1]):
                if nbr[j, kk] == i:
                    rev[i, k] = kk
                    break
    return rev


@njit(cache=True)
def edge_kappa(th, domain, L_km, nbr, dist, axis, tab, T_ax, L_ax, S_ax, airy, k_normal, k_max, pinned=None,
               tab_pinned=None, rev=None):
    """Per-edge diffusivity (m^2/yr): kappa_normal + min(q / |grad T|, k_max), q from the table at the edge's mean T and
    L and its slope (from tab_pinned on edges with a pinned node, when given). With rev (reverse_index), each edge is
    evaluated once (from its lower node) and written to both of its slots, so K is exactly symmetric. Returns (K (n, 8),
    largest K, edges at k_max, edges counted)."""
    n = th.shape[0]
    K = np.zeros(nbr.shape)
    kmax = 0.
    capped = 0
    counted = 0
    smin = 10.**S_ax[0]
    for i in range(n):
        if not domain[i]:
            continue
        for kk in range(nbr.shape[1]):
            j = nbr[i, kk]
            if j < 0 or not axis[i, kk] or not domain[j]:
                continue
            if rev is not None and (j < i or rev[i, kk] < 0):
                continue
            tb = tab
            if pinned is not None and (pinned[i] or pinned[j]):
                tb = tab_pinned
            kf = _edge_k(th[i], th[j], dist[i, kk], L_km[i], L_km[j], tb, T_ax, L_ax, S_ax, airy, smin)
            counted += 1
            if kf >= k_max:
                kf = k_max
                capped += 1
            K[i, kk] = k_normal+kf
            if rev is not None:
                K[j, rev[i, kk]] = K[i, kk]
            if K[i, kk] > kmax:
                kmax = K[i, kk]
    return K, kmax, capped, counted


@njit(cache=True)
def diffuse_edges(z, domain, K, dt, nbr, dist, axis, cell_area):
    """graph_landscape.diffuse_variable with a per-edge diffusivity (symmetric: K[i, k] = K[j, k'] for the same edge,
    since it is computed from the edge's symmetric mean and slope); conservative."""
    n = z.shape[0]
    out = z.copy()
    for i in range(n):
        if not domain[i]:
            continue
        s = 0.
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0 or not axis[i, k] or not domain[j]:
                continue
            d = dist[i, k]
            s += K[i, k]*(z[j]-z[i])/d*(np.sqrt(cell_area[i]*cell_area[j])/d)
        out[i] = z[i]+dt*s/cell_area[i]
    return out


@njit(cache=True)
def hot_subcycle(th, hot, hot_nodes, L_km, nbr, dist, axis, cell_area, tab, T_ax, L_ax, S_ax, airy, k_normal, k_max,
                 pinned, tab_pinned, use_pinned, h_total, bound, rev):
    """Integrate the hot edges (hot[i, k]: the edge's diffusivity above the split) over h_total with their own explicit
    substeps, re-evaluating only their diffusivity; the other edges are left to the caller's step (operator split,
    each part stable and conservative). Returns (th, substeps)."""
    smin = 10.**S_ax[0]
    K = np.zeros(nbr.shape)
    rem = h_total
    steps = 0
    out = th.copy()
    while rem > 0.:
        kmax = 0.
        for ii in range(hot_nodes.shape[0]):
            i = hot_nodes[ii]
            for kk in range(nbr.shape[1]):
                if not hot[i, kk]:
                    continue
                j = nbr[i, kk]
                if j < i:
                    continue
                tb = tab
                if use_pinned and (pinned[i] or pinned[j]):
                    tb = tab_pinned
                kf = _edge_k(out[i], out[j], dist[i, kk], L_km[i], L_km[j], tb, T_ax, L_ax, S_ax, airy, smin)
                if kf > k_max:
                    kf = k_max
                K[i, kk] = k_normal+kf
                K[j, rev[i, kk]] = K[i, kk]
                if K[i, kk] > kmax:
                    kmax = K[i, kk]
        if kmax <= 0.:
            break
        h = min(rem, .9/(kmax*bound))
        new = np.empty(hot_nodes.shape[0])
        for ii in range(hot_nodes.shape[0]):
            i = hot_nodes[ii]
            s = 0.
            for kk in range(nbr.shape[1]):
                if not hot[i, kk]:
                    continue
                j = nbr[i, kk]
                d = dist[i, kk]
                s += K[i, kk]*(out[j]-out[i])/d*(np.sqrt(cell_area[i]*cell_area[j])/d)
            new[ii] = out[i]+h*s/cell_area[i]
        for ii in range(hot_nodes.shape[0]):
            out[hot_nodes[ii]] = new[ii]
        rem -= h
        steps += 1
    return out, steps


def lith_km(sim):
    """L (km) per node from C's thermal age; the frame's reference where C is off or the age is unknown."""
    from planet_tectonics.lithosphere import frame, thickness_km
    p, w = sim.p, sim.world
    ref = frame(p)[0]
    if p['lithosphere'] != 'thermal_age':
        return np.full(sim.g.size, float(ref))
    age = sim.t-w['therm']
    return np.where(np.isfinite(age), thickness_km(np.nan_to_num(age, nan=0.), p), ref)


def pinned_nodes(sim):
    """'rey': nodes within flow_pinned_km of active convergence (the arc belt over a converging trench; collision
    boundary nodes closing faster than boundary_min_mm_yr), from the last step's boundaries."""
    from planet_landscape import graph_landscape as gl
    from planet_tectonics.params import COLLISION
    g, p = sim.g, sim.p
    active = np.asarray(getattr(sim, 'belt', np.zeros(g.size, bool)), bool).copy()
    if hasattr(sim, 'btype') and hasattr(sim, 'conv_rate'):
        active |= (sim.btype == COLLISION) & (sim.conv_rate*1e3 > p['boundary_min_mm_yr'])
    if not active.any():
        return active
    reach = p['flow_pinned_km']*1e3
    return gl.distance_to(active, g.nbr, g.dist, reach) <= reach


def flow(sim, th, dt_myr):
    """One step of the crust flow with the measured rheology on continental nodes. Returns the new thickness and books
    the substeps and the share of edges at the numerical ceiling in sim.erosion_totals."""
    g, p, w = sim.g, sim.p, sim.world
    pinned = tab_p = None
    if p['crust_flow_rule'] == 'rey':
        tab, tab_p = flux_table(p, 'sheet'), flux_table(p, 'channel')
        pinned = pinned_nodes(sim)
        sim.erosion_totals['crust_flow_pinned_area'] = sim.erosion_totals.get('crust_flow_pinned_area', 0.)+float(
            g.area[pinned & w['cont']].sum()/g.area[w['cont']].sum())*dt_myr
    else:
        tab = flux_table(p)
    L = lith_km(sim)
    airy = 1.-p['rho_crust']/p['rho_mantle']
    k_split = p['flow_kappa_split_m2_yr']
    use_p = pinned is not None
    pin = pinned if use_p else np.zeros(1, bool)
    tp = tab_p if use_p else tab
    T = dt_myr*1e6
    tot = sim.erosion_totals
    if getattr(sim, '_flow_rev', None) is None:
        sim._flow_rev = reverse_index(g.nbr)
    rev = sim._flow_rev
    while T > 0:
        K, kmax, capped, counted = edge_kappa(th, w['cont'], L, g.nbr, g.dist, g.axis, tab, T_AXIS_KM, L_AXIS_KM,
                                              LOGS_AXIS, airy, p['crust_kappa_normal_m2_yr'], p['flow_kappa_max_m2_yr'],
                                              pinned, tab_p, rev)
        if kmax <= 0:
            break
        # operator split (numerical): edges above k_split take their own substeps; the rest one step of h
        hot = K > k_split
        h = min(T, .9/(min(kmax, k_split)*sim.bound), p['flow_split_step_myr']*1e6)
        if hot.any():
            hn = np.flatnonzero(hot.any(1))
            th, inner = hot_subcycle(th, hot, hn, L, g.nbr, g.dist, g.axis, g.area, tab, T_AXIS_KM, L_AXIS_KM, LOGS_AXIS,
                                     airy, p['crust_kappa_normal_m2_yr'], p['flow_kappa_max_m2_yr'], pin, tp, use_p, h,
                                     sim.bound, rev)
            tot['crust_flow_substeps'] = tot.get('crust_flow_substeps', 0.)+float(inner)
            K = np.where(hot, 0., K)
        th = diffuse_edges(th, w['cont'], K, h, g.nbr, g.dist, g.axis, g.area)
        tot['crust_flow_substeps'] = tot.get('crust_flow_substeps', 0.)+1.
        tot['crust_flow_capped_edges'] = tot.get('crust_flow_capped_edges', 0.)+float(capped)*h/1e6
        tot['crust_flow_edges'] = tot.get('crust_flow_edges', 0.)+float(counted)*h/1e6
        T -= h
    return th
