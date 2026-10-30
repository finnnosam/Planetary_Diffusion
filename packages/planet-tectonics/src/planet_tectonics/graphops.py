"""Graph kernels on the cube-sphere node graph: nearest labelled node, connected components, holes in a plate mask."""
import heapq
import numpy as np

from planet_grid.jit import njit



@njit(cache=True)
def nearest_label(source_label, nbr, dist):
    """Label of the nearest labelled node (label >= 0) for every node (multi-source Dijkstra)."""
    n = source_label.shape[0]
    d = np.full(n, 1e30)
    lab = source_label.copy()
    heap = [(0.0, -1)]
    heap.pop()
    for i in range(n):
        if source_label[i] >= 0:
            d[i] = 0.
            heapq.heappush(heap, (0., i))
    while len(heap) > 0:
        h, i = heapq.heappop(heap)
        if h > d[i]:
            continue
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0:
                continue
            v = h+dist[i, k]
            if v < d[j]:
                d[j] = v; lab[j] = lab[i]
                heapq.heappush(heap, (v, j))
    return lab, d



@njit(cache=True)
def components(mask, nbr):
    """Connected-component label of every mask node (-1 elsewhere)."""
    n = mask.shape[0]
    lab = np.full(n, -1, np.int64)
    work = np.empty(n, np.int64)
    c = 0
    for s in range(n):
        if not mask[s] or lab[s] >= 0:
            continue
        lab[s] = c; work[0] = s; w = 1
        while w > 0:
            w -= 1; i = work[w]
            for k in range(nbr.shape[1]):
                j = nbr[i, k]
                if j >= 0 and mask[j] and lab[j] < 0:
                    lab[j] = c; work[w] = j; w += 1
        c += 1
    return lab



@njit(cache=True)
def plate_holes(mask, nbr):
    """Unset raster nodes whose every neighbour is set, and the neighbour to copy from (sampling holes left by the
    nearest-node forward/backward maps not being exact inverses; real plate edges always touch unset nodes)."""
    n = mask.shape[0]
    hole = np.zeros(n, np.bool_)
    src = np.full(n, -1, np.int64)
    for i in range(n):
        if mask[i]:
            continue
        ok = True; first = -1
        for k in range(nbr.shape[1]):
            j = nbr[i, k]
            if j < 0:
                continue
            if not mask[j]:
                ok = False; break
            if first < 0:
                first = j
        if ok and first >= 0:
            hole[i] = True; src[i] = first
    return hole, src
