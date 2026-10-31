"""Boundary typing: ridge, subduction, collision or transform at each boundary node, from the plates' relative motion
(node rule and two-ring rule)."""
import numpy as np

from planet_grid.jit import njit



@njit(cache=True)
def classify_boundaries(owner, slot, omega, xyz, nbr, score, cont, radius, vmin):
    """Boundary type, convergence rate (m/yr, + = closing) and trench flag of every node. A boundary whose normal
    (closing or opening) rate is under half its shear rate, or under vmin, is a transform.

    omega: (K, 3) angular velocity vectors (rad/yr) by plate slot; slot[owner] gives the row.
    A node is a trench node when it would subduct under its converging neighbour (lower score)."""
    n = owner.shape[0]
    btype = np.zeros(n, np.int8)
    rate = np.zeros(n)
    trench = np.zeros(n, np.bool_)
    for i in range(n):
        a = slot[owner[i]]
        best = 0.
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0 or owner[j] == owner[i]:
                continue
            b = slot[owner[j]]
            w0 = omega[a, 0]-omega[b, 0]; w1 = omega[a, 1]-omega[b, 1]; w2 = omega[a, 2]-omega[b, 2]
            x0, x1, x2 = xyz[i, 0], xyz[i, 1], xyz[i, 2]
            v0 = (w1*x2-w2*x1)*radius; v1 = (w2*x0-w0*x2)*radius; v2 = (w0*x1-w1*x0)*radius
            d0 = xyz[j, 0]-x0; d1 = xyz[j, 1]-x1; d2 = xyz[j, 2]-x2
            dn = d0*x0+d1*x1+d2*x2
            d0 -= dn*x0; d1 -= dn*x1; d2 -= dn*x2
            L = np.sqrt(d0*d0+d1*d1+d2*d2)
            conv = (v0*d0+v1*d1+v2*d2)/L
            vt = np.sqrt(max(v0*v0+v1*v1+v2*v2-conv*conv, 0.))   # boundary-parallel (shear) component
            if abs(conv) >= abs(best) or btype[i] == 0:
                best = conv
                if abs(conv) < max(vmin, .5*vt):
                    btype[i] = 4; trench[i] = False
                elif conv > 0:
                    if cont[i] and cont[j]:
                        btype[i] = 3
                        trench[i] = False
                    else:
                        btype[i] = 2
                        trench[i] = score[i] < score[j]
                else:
                    btype[i] = 1; trench[i] = False
        rate[i] = best
    return btype, rate, trench



@njit(cache=True)
def classify_boundaries_ring(owner, slot, omega, xyz, nbr, score, cont, radius, vmin):
    """As classify_boundaries, with the boundary normal from the node's two-ring (PARAMS 'boundary_normal' 'ring'):
    toward each neighbouring plate b, the tangent direction of sum_k (x_k - x_i) over the node's neighbours and their
    neighbours, +1 on b and -1 on the node's own plate (other plates 0). Of the neighbouring plates, the one with the
    largest |normal rate| types the node, as before; a trench node is the one with the lower score."""
    n = owner.shape[0]
    K = nbr.shape[1]
    btype = np.zeros(n, np.int8)
    rate = np.zeros(n)
    trench = np.zeros(n, np.bool_)
    ring = np.empty(K*K+K, np.int64)
    for i in range(n):
        a = owner[i]
        x0, x1, x2 = xyz[i, 0], xyz[i, 1], xyz[i, 2]
        nr = -1
        best = 0.
        done = False
        for k in range(K):
            j = nbr[i, k]
            if j < 0 or owner[j] == a:
                continue
            b = owner[j]
            seen = False                                   # each neighbouring plate once
            for k0 in range(k):
                j0 = nbr[i, k0]
                if j0 >= 0 and owner[j0] == b:
                    seen = True
            if seen:
                continue
            if nr < 0:                                     # the two-ring, built once per boundary node
                nr = 0
                for k1 in range(K):
                    j1 = nbr[i, k1]
                    if j1 < 0:
                        continue
                    ring[nr] = j1; nr += 1
                    for k2 in range(K):
                        j2 = nbr[j1, k2]
                        if j2 >= 0 and j2 != i:
                            ring[nr] = j2; nr += 1
            d0 = 0.; d1 = 0.; d2 = 0.
            for r in range(nr):
                q = ring[r]
                sg = 1. if owner[q] == b else (-1. if owner[q] == a else 0.)
                d0 += sg*(xyz[q, 0]-x0); d1 += sg*(xyz[q, 1]-x1); d2 += sg*(xyz[q, 2]-x2)
            dn = d0*x0+d1*x1+d2*x2
            d0 -= dn*x0; d1 -= dn*x1; d2 -= dn*x2
            L = np.sqrt(d0*d0+d1*d1+d2*d2)
            if L == 0.:
                continue
            sa, sb = slot[a], slot[b]
            w0 = omega[sa, 0]-omega[sb, 0]; w1 = omega[sa, 1]-omega[sb, 1]; w2 = omega[sa, 2]-omega[sb, 2]
            v0 = (w1*x2-w2*x1)*radius; v1 = (w2*x0-w0*x2)*radius; v2 = (w0*x1-w1*x0)*radius
            conv = (v0*d0+v1*d1+v2*d2)/L
            vt = np.sqrt(max(v0*v0+v1*v1+v2*v2-conv*conv, 0.))   # boundary-parallel (shear) component
            if abs(conv) >= abs(best) or not done:
                best = conv
                done = True
                if abs(conv) < max(vmin, .5*vt):
                    btype[i] = 4; trench[i] = False
                elif conv > 0:
                    if cont[i] and cont[j]:
                        btype[i] = 3
                        trench[i] = False
                    else:
                        btype[i] = 2
                        trench[i] = score[i] < score[j]
                else:
                    btype[i] = 1; trench[i] = False
        rate[i] = best
    return btype, rate, trench



@njit(cache=True)
def boundary_frame_ring(owner, slot, omega, xyz, nbr, radius):
    """For every boundary node, the neighbouring plate that types it under classify_boundaries_ring (the one with the
    largest |normal rate|) and the unit boundary normal toward that plate (the two-ring normal, tangent at the node);
    -1 and zeros elsewhere. Used by force_law 'torque' (forty-third pass) for the direction of boundary forces."""
    n = owner.shape[0]
    K = nbr.shape[1]
    other = np.full(n, -1, np.int64)
    normal = np.zeros((n, 3))
    ring = np.empty(K*K+K, np.int64)
    for i in range(n):
        a = owner[i]
        x0, x1, x2 = xyz[i, 0], xyz[i, 1], xyz[i, 2]
        nr = -1
        best = 0.
        done = False
        for k in range(K):
            j = nbr[i, k]
            if j < 0 or owner[j] == a:
                continue
            b = owner[j]
            seen = False
            for k0 in range(k):
                j0 = nbr[i, k0]
                if j0 >= 0 and owner[j0] == b:
                    seen = True
            if seen:
                continue
            if nr < 0:
                nr = 0
                for k1 in range(K):
                    j1 = nbr[i, k1]
                    if j1 < 0:
                        continue
                    ring[nr] = j1; nr += 1
                    for k2 in range(K):
                        j2 = nbr[j1, k2]
                        if j2 >= 0 and j2 != i:
                            ring[nr] = j2; nr += 1
            d0 = 0.; d1 = 0.; d2 = 0.
            for r in range(nr):
                q = ring[r]
                sg = 1. if owner[q] == b else (-1. if owner[q] == a else 0.)
                d0 += sg*(xyz[q, 0]-x0); d1 += sg*(xyz[q, 1]-x1); d2 += sg*(xyz[q, 2]-x2)
            dn = d0*x0+d1*x1+d2*x2
            d0 -= dn*x0; d1 -= dn*x1; d2 -= dn*x2
            L = np.sqrt(d0*d0+d1*d1+d2*d2)
            if L == 0.:
                continue
            sa, sb = slot[a], slot[b]
            w0 = omega[sa, 0]-omega[sb, 0]; w1 = omega[sa, 1]-omega[sb, 1]; w2 = omega[sa, 2]-omega[sb, 2]
            v0 = (w1*x2-w2*x1)*radius; v1 = (w2*x0-w0*x2)*radius; v2 = (w0*x1-w1*x0)*radius
            conv = (v0*d0+v1*d1+v2*d2)/L
            if abs(conv) >= abs(best) or not done:
                best = conv
                done = True
                other[i] = b
                normal[i, 0] = d0/L; normal[i, 1] = d1/L; normal[i, 2] = d2/L
    return other, normal
