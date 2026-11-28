"""Small helpers shared by the planet layer: node locations, random spherical fields, rotations, fast cube lookup."""
from functools import cache

import numpy as np

from planet_grid import cube
from planet_grid.jit import njit, prange


def node_loc(xyz):
    """(N, 3) unit vectors -> (N, 2) colatitude, longitude in [0, 2 pi) (ducc0 convention)."""
    x, y, z = np.asarray(xyz, np.float64).T
    return np.ascontiguousarray(np.stack((np.arccos(np.clip(z, -1, 1)), np.mod(np.arctan2(y, x), 2*np.pi)), 1))


def random_field(xyz, lmax, slope, rng):
    """Zero-mean, unit-variance random field with angular power ~ l^-slope (1 <= l <= lmax), at the given points."""
    import ducc0
    alm = np.zeros((lmax+1)*(lmax+2)//2, np.complex128)
    for m in range(lmax+1):
        for l in range(max(m, 1), lmax+1):
            i = m*(2*lmax+1-m)//2+l
            a = rng.normal()*l**(-slope/2)
            alm[i] = a if m == 0 else (a+1j*rng.normal()*l**(-slope/2))/np.sqrt(2)
    f = ducc0.sht.experimental.synthesis_general(alm=alm[None], spin=0, lmax=lmax, loc=node_loc(xyz), epsilon=1e-8)[0]
    return (f-f.mean())/f.std()


def spectrum_field(xyz, power, rng):
    """Gaussian isotropic random field with a given expected power per degree, at the given points, not rescaled.

    power[l] (l = 0 .. lmax; power[0] is ignored, so the field has no degree-0 part) is the expected sum over m of
    the squared coefficients on orthonormal real harmonics (integral of Y^2 over the sphere = 1), in the field's
    units squared: each of the 2l + 1 coefficients has variance power[l] / (2l + 1), and the expected mean square
    over the sphere is sum(power[1:]) / (4 pi). One realization's power scatters around that (chi-square with 2l + 1
    degrees of freedom per degree)."""
    import ducc0
    power = np.asarray(power, np.float64)
    lmax = power.size-1
    alm = np.zeros((lmax+1)*(lmax+2)//2, np.complex128)
    for m in range(lmax+1):
        for l in range(max(m, 1), lmax+1):
            i = m*(2*lmax+1-m)//2+l
            s = np.sqrt(power[l]/(2*l+1))
            # ducc0's real-field convention: f = sum_l a_l0 Y_l0 + 2 Re sum_{m>0} a_lm Y_lm, so a real orthonormal
            # coefficient of variance s^2 is a_lm = (x + i y) / sqrt(2) with x, y ~ N(0, s^2)
            alm[i] = rng.normal()*s if m == 0 else (rng.normal()+1j*rng.normal())*s/np.sqrt(2)
    return ducc0.sht.experimental.synthesis_general(alm=alm[None], spin=0, lmax=lmax, loc=node_loc(xyz), epsilon=1e-8)[0]


def rotation(axis, angle):
    """Rotation matrix about a unit axis (right-handed, radians)."""
    a = np.asarray(axis, np.float64); a = a/np.linalg.norm(a)
    k = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3)+np.sin(angle)*k+(1-np.cos(angle))*(k@k)


_NORMAL, _RIGHT, _DOWN = cube.NORMAL.copy(), cube.RIGHT.copy(), cube.DOWN.copy()
_FACE_POS = np.array([0, 1, 4]); _FACE_NEG = np.array([2, 3, 5])
_CHART_BINS = 16384


@cache
def _chart_table(n):
    """A rounded chart index per ratio bin, with -1 wherever rounding can change inside the bin."""
    edges = np.linspace(-1., 1., _CHART_BINS+1)
    # Pad for the division and bin scaling, then for library atan rounding and the chart arithmetic.
    eps = np.finfo(np.float64).eps
    lower = n*(1+np.arctan(edges[:-1]-16*eps)/cube.QUARTER)/2 - 32*n*eps
    upper = n*(1+np.arctan(edges[1:]+16*eps)/cube.QUARTER)/2 + 32*n*eps
    lo = np.floor(lower+.5)
    hi = np.floor(upper+.5)
    table = np.where(lo == hi, lo, -1).astype(np.int32)
    table.flags.writeable = False
    return table


@njit(cache=True, inline='always')
def _chart_index(num, den, n, table):
    if table is not None and 0. < den < np.inf:
        # The selected face has the largest component, so num/den is in [-1, 1].
        bin_index = min(int((num/den+1.)*(_CHART_BINS//2)), _CHART_BINS-1)
        index = table[bin_index]
        if index >= 0:
            return index
    x = n*(1+np.arctan2(num, den)/(np.pi/4))/2
    return min(max(int(np.floor(x+.5)), 0), n)


@njit(cache=True, parallel=True)
def nearest_entries(p, n, normal, right, down, face_pos, face_neg, table=None):
    """Shared-layout cube entry of the node nearest (in chart index) to each direction p (M, 3); as cube.coordinates.
    Each point on its own (prange: bit-identical at any thread count)."""
    m = p.shape[0]
    out = np.empty(m, np.int64)
    for i in prange(m):
        ax = 0
        if abs(p[i, 1]) > abs(p[i, ax]):
            ax = 1
        if abs(p[i, 2]) > abs(p[i, ax]):
            ax = 2
        f = face_pos[ax] if p[i, ax] >= 0 else face_neg[ax]
        den = p[i, 0]*normal[f, 0]+p[i, 1]*normal[f, 1]+p[i, 2]*normal[f, 2]
        r = p[i, 0]*right[f, 0]+p[i, 1]*right[f, 1]+p[i, 2]*right[f, 2]
        d = p[i, 0]*down[f, 0]+p[i, 1]*down[f, 1]+p[i, 2]*down[f, 2]
        xi = _chart_index(r, den, n, table)
        yi = _chart_index(d, den, n, table)
        out[i] = (f*(n+1)+yi)*(n+1)+xi
    return out


def nearest_nodes(graph, p):
    """Graph node nearest to each direction p (M, 3), for a planet_grid.sphere_graph.CubeSphere."""
    n = graph.n
    e = nearest_entries(np.ascontiguousarray(p, np.float64), n, _NORMAL, _RIGHT, _DOWN, _FACE_POS, _FACE_NEG,
                        _chart_table(n))
    return graph.node_of_entry[e]


def equirect_nodes(graph, height):
    """Node index for every pixel of a (height, 2*height) equirectangular grid (nearest node; categorical-safe)."""
    lat = np.pi/2-np.pi*(np.arange(height)+.5)/height
    lon = -np.pi+2*np.pi*(np.arange(2*height)+.5)/(2*height)
    la, lo = np.meshgrid(lat, lon, indexing='ij')
    p = np.stack((np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)), -1).reshape(-1, 3)
    return nearest_nodes(graph, p).reshape(height, 2*height)
