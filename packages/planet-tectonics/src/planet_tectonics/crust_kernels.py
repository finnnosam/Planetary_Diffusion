"""Crust kernels: start-crust quantile mapping, margin shortening, edge collapse, pure-shear rifting, arc-band spreading
and trench fluxes from consumed area."""
import heapq
import numpy as np

from planet_grid.jit import njit

from planet_tectonics.graphops import nearest_label
from planet_tectonics.params import EARTH_STABLE_CRUST_KM, EARTH_STABLE_CRUST_Q



def quantile_map(values, weights, q=EARTH_STABLE_CRUST_Q, target=EARTH_STABLE_CRUST_KM):
    """Map values onto a target distribution by their area-weighted rank (mid-rank of each node's area), keeping their
    order; ranks beyond the outermost quantiles take the end values."""
    o = np.argsort(values, kind='stable')
    c = np.cumsum(weights[o]); u = (c-.5*weights[o])/c[-1]
    out = np.empty(len(values))
    out[o] = np.interp(u, q, target)
    return out



def margin_shortening(near, d, dom, s, thick, area, dt, d_belt, sg_belt, d_fore, sg_fore, cap=.3):
    """Volume-conserving upper-plate shortening (m of crust per node over dt). For every trench segment (the nodes
    whose nearest overriding trench node is the same, `near`), crust is taken from the foreland band (Gaussian in
    trench distance d at d_fore +- sg_fore) at s * thick per unit trench length, at most `cap` of a column, and the
    same volume is added to the belt (Gaussian at d_belt +- sg_belt) of that segment. dom: continental nodes of the
    overriding plate within reach. Segments without belt nodes move nothing."""
    gf = np.exp(-.5*((d-d_fore)/sg_fore)**2)/(sg_fore*np.sqrt(2*np.pi))
    gb = np.exp(-.5*((d-d_belt)/sg_belt)**2)/(sg_belt*np.sqrt(2*np.pi))
    seg = np.where(dom, near, -1)
    m = seg >= 0
    nseg = int(near.max())+1 if m.any() else 0
    out = np.zeros(d.shape[0])
    if nseg == 0:
        return out
    bw = np.where(m, gb*area, 0.)
    B = np.bincount(seg[m], weights=bw[m], minlength=nseg)
    sink = np.where(m, -np.minimum(s*dt*thick*gf, cap*thick), 0.)
    sink = np.where(m & (B[np.maximum(seg, 0)] > 0), sink, 0.)
    Rv = np.bincount(seg[m], weights=(sink*area)[m], minlength=nseg)            # m^3 removed per segment (<= 0)
    add = np.where(m, -Rv[np.maximum(seg, 0)]*bw/np.maximum(B[np.maximum(seg, 0)], 1e-300)/area, 0.)
    return sink+add



@njit(cache=True)
def edge_spread(th, cont, kap, excess, owner, ocean_ok, h, nbr, dist, axis, area):
    """Collapse flux of thick continental crust across the coastline (one explicit substep of length h, years): from a
    continental node i with excess thickness E_i over the collapse threshold into each axis neighbour j that is
    oceanic, owned by the same plate and not on a plate boundary: kap_i * E_i / d * (sqrt(a_i a_j) / d) * h m^3.
    Conservative; never takes more than the excess."""
    n = th.shape[0]
    out = th.copy()
    for i in range(n):
        if not cont[i] or excess[i] <= 0. or kap[i] <= 0.:
            continue
        tot = 0.
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0 or not axis[i, k] or cont[j] or not ocean_ok[j] or owner[j] != owner[i]:
                continue
            d = dist[i, k]
            tot += kap[i]*excess[i]/d*(np.sqrt(area[i]*area[j])/d)*h
        if tot <= 0.:
            continue
        scale = min(1., excess[i]*area[i]/tot)
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0 or not axis[i, k] or cont[j] or not ocean_ok[j] or owner[j] != owner[i]:
                continue
            d = dist[i, k]
            f = kap[i]*excess[i]/d*(np.sqrt(area[i]*area[j])/d)*h*scale
            out[i] -= f/area[i]
            out[j] += f/area[j]
    return out



@njit(cache=True)
def _zone_dijkstra(g, zone, nbr, dist, cutoff, dd, touched):
    """Distances from node g to the zone nodes within cutoff, through zone nodes only. Fills dd (inf elsewhere) and
    returns the number of touched nodes listed in `touched` (g first); the caller resets dd[touched[:m]] = inf."""
    heap = [(0.0, g)]
    dd[g] = 0.
    touched[0] = g; m = 1
    while len(heap) > 0:
        h, i = heapq.heappop(heap)
        if h > dd[i]:
            continue
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0 or not zone[j]:
                continue
            v = h+dist[i, k]
            if v < dd[j] and v <= cutoff:
                if dd[j] == np.inf:
                    touched[m] = j; m += 1
                dd[j] = v
                heapq.heappush(heap, (v, j))
    return m



@njit(cache=True)
def rift_pure_shear(rows, area, thick, zone, nbr, dist, W, Tb):
    """Pure-shear stretching of the zone around rift gap nodes (rift_rule 'pure_shear', see PARAMS).
    rows: gap nodes; zone: continental nodes that may stretch (not gaps); W: half-width (m) of the linear strain taper,
    measured from the gap cell's edge (half its spacing). Returns per row the strain-weighted mean zone thickness Tbar
    (before stretching), the gap crust Tg (0 where it would be thinner than Tb: breakup, or no zone), and per node the
    thinning dT (<= 0).
    Zone node i takes c_i = sum over rifting gaps g of a_g k_ig / K_g (K_g = sum_j a_j k_jg) and thins to
    T_i exp(-c_i); gap g receives its share a_g k_ig / (K_g c_i) of every removed volume, so the volume is conserved
    exactly; Tg = Tbar for small strain, less for a gap that is wide against W (a 46 km cell opening into a 100 km
    zone is a finite stretch)."""
    n = area.shape[0]
    R = rows.shape[0]
    dd = np.full(n, np.inf)
    touched = np.empty(n, np.int64)
    K = np.zeros(R); S = np.zeros(R)
    for r in range(R):
        g = rows[r]
        rg = .5*np.sqrt(area[g])
        m = _zone_dijkstra(g, zone, nbr, dist, W+rg, dd, touched)
        for q in range(1, m):
            j = touched[q]
            k = 1.-max(dd[j]-rg, 0.)/W
            if k > 0.:
                K[r] += area[j]*k; S[r] += area[j]*k*thick[j]
        for q in range(m):
            dd[touched[q]] = np.inf
    Tbar = np.where(K > 0., S/np.maximum(K, 1e-300), 0.)
    ok = (K > 0.) & (Tbar >= Tb)
    # finite strain makes Tg < Tbar (the gap is part of the stretched zone); gaps whose Tg falls below Tb break up,
    # and the rest is recomputed without them (a few rounds at most)
    for it in range(20):
        c, dT, Tg = _pure_shear_apply(rows, ok, K, area, thick, zone, nbr, dist, W, dd, touched)
        bad = ok & (Tg < Tb)
        if not bad.any():
            break
        ok &= ~bad
    return Tbar, Tg, dT



@njit(cache=True)
def _pure_shear_apply(rows, ok, K, area, thick, zone, nbr, dist, W, dd, touched):
    n = area.shape[0]
    R = rows.shape[0]
    c = np.zeros(n)
    for r in range(R):
        if not ok[r]:
            continue
        g = rows[r]
        rg = .5*np.sqrt(area[g])
        m = _zone_dijkstra(g, zone, nbr, dist, W+rg, dd, touched)
        for q in range(1, m):
            j = touched[q]
            k = 1.-max(dd[j]-rg, 0.)/W
            if k > 0.:
                c[j] += area[g]*k/K[r]
        for q in range(m):
            dd[touched[q]] = np.inf
    dT = np.zeros(n)
    for i in range(n):
        if c[i] > 0.:
            dT[i] = -thick[i]*(1.-np.exp(-c[i]))
    Tg = np.zeros(R)
    for r in range(R):
        if not ok[r]:
            continue
        g = rows[r]
        rg = .5*np.sqrt(area[g])
        m = _zone_dijkstra(g, zone, nbr, dist, W+rg, dd, touched)
        v = 0.
        for q in range(1, m):
            j = touched[q]
            k = 1.-max(dd[j]-rg, 0.)/W
            if k > 0.:
                v += -dT[j]*area[j]*(area[g]*k/K[r])/c[j]
        for q in range(m):
            dd[touched[q]] = np.inf
        Tg[r] = v/area[g]
    return c, dT, Tg



def spread_to_band(vol, near, d, band, d0, sg, area):
    """Thickness (m) that spreads the volumes vol (m^3 at nodes) over the arc band of each one's nearest trench node:
    Gaussian in trench distance d at d0 +- sg over the band nodes with the same nearest trench node, normalized so the
    volume is conserved. A volume whose trench node has no band node stays where it is (B-d, thirty-ninth pass)."""
    n = vol.size
    out = np.zeros(n)
    src = np.flatnonzero(vol > 0)
    if src.size == 0:
        return out
    key = np.where(near >= 0, near, -1)
    V = np.zeros(n); np.add.at(V, key[src][key[src] >= 0], vol[src][key[src] >= 0])
    wgt = np.where(band & (key >= 0), np.exp(-.5*((d-d0)/sg)**2), 0.)
    S = np.zeros(n); np.add.at(S, np.where(wgt > 0, key, 0), wgt*area*(wgt > 0))
    has = (wgt > 0) & (S[np.clip(key, 0, None)] > 0)
    out[has] = V[key[has]]*wgt[has]/S[key[has]]
    lone = (key[src] < 0) | (S[np.clip(key[src], 0, None)] <= 0)          # no band: stays on its node
    out[src[lone]] += vol[src[lone]]/area[src[lone]]
    return out



def trench_flux_consumed(A, near, d, ok, terms, area, min_len, own, nbr, dist):
    """Trench fluxes from consumed area (trench_flux 'consumed', B-d2, fortieth pass). A: oceanic area (m^2) consumed
    at each world node last step; near / d: nearest overriding trench node and the distance to it; ok: nodes of that
    node's plate within reach. Each consumed area is routed to the nearest overriding trench node of its plate that owns
    a band (Gaussian integral at arc_d0 150 +- 50 km over its band >= min_len; on a staircase boundary the outer row owns
    almost none), only where ok (consumption with no overriding trench of the receiving plate in reach, e.g. overlap at a
    transform, feeds nothing). own / nbr / dist: plate of each node and the graph. terms: (mask, H, d0, sg); each
    adds A_s x H on its mask's nodes of the band of s, weighted by its Gaussian in d and normalized by the Gaussian's
    integral over the whole band of s (all of its ok nodes), floored at min_len (m; a truncated band keeps the share it
    has, as the 'band' law does). Returns the thickness change (m) and the routed and unrouted areas (m^2)."""
    n = area.size
    out = np.zeros(n)
    if A is None or not (A > 0).any():
        return out, 0., 0.
    key = np.where(ok, near, -1)
    ref = np.exp(-.5*((d-150e3)/50e3)**2)/(50e3*np.sqrt(2*np.pi))
    owns_band = np.bincount(key[ok], weights=(ref*area)[ok], minlength=n) >= min_len
    near2, _ = nearest_label(np.where(owns_band, np.arange(n), -1).astype(np.int64), nbr, dist)
    to = np.where(ok & (near2 >= 0), near2, -1)
    to[to >= 0] = np.where(own[to[to >= 0]] == own[to >= 0], to[to >= 0], -1)
    got = (A > 0) & (to >= 0)
    As = np.bincount(to[got], weights=A[got], minlength=n)
    routed, unrouted = float(A[got].sum()), float(A[(A > 0) & ~got].sum())
    kk = np.clip(key, 0, None)
    for mask, H, d0, sg in terms:
        if H == 0:
            continue
        shape = np.exp(-.5*((d-d0)/sg)**2)/(sg*np.sqrt(2*np.pi))
        norm = np.bincount(key[ok], weights=(shape*area)[ok], minlength=n)
        norm = np.maximum(norm, min_len)
        out = out+np.where(mask & ok, As[kk]*H*shape/norm[kk], 0.)
    return out, routed, unrouted
