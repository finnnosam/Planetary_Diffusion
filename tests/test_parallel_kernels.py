"""The parallel kernels and the step caches leave results bit-identical (journal 154)."""
import numpy as np

from planet_grid import cube, fields
from planet_grid.jit import set_threads
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import crust_flow
from planet_tectonics.params import PARAMS
from planet_tectonics.plates import PlateSim


def _at_threads(fn, threads):
    out = []
    for t in threads:
        set_threads(t)
        out.append(fn())
    set_threads(None)
    return out


def test_nearest_nodes_and_edge_kappa_do_not_depend_on_threads():
    g = CubeSphere(64)
    rng = np.random.default_rng(3)
    p = rng.normal(size=(50000, 3))
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    a, b = _at_threads(lambda: fields.nearest_nodes(g, p), (1, 4))
    assert np.array_equal(a, b)
    th = 30e3+20e3*rng.random(g.size)
    dom = rng.random(g.size) < .6
    L = 100.+100.*rng.random(g.size)
    tab = crust_flow.flux_table({**PARAMS, 'crust_flow_rule': 'sheet'})
    rev = crust_flow.reverse_index(g.nbr)
    args = (th, dom, L, g.nbr, g.dist, g.axis, tab, crust_flow.T_AXIS_KM, crust_flow.L_AXIS_KM, crust_flow.LOGS_AXIS,
            .18, 1e2, 1e7)
    r1, r4 = _at_threads(lambda: crust_flow.edge_kappa(*args, None, None, rev), (1, 4))
    assert np.array_equal(r1[0], r4[0]) and r1[1:] == r4[1:]
    ok = rev >= 0                                      # each edge written once, to both of its slots
    assert np.array_equal(r1[0][ok], r1[0][g.nbr[ok], rev[ok]])
    d1, d4 = _at_threads(lambda: crust_flow.diffuse_edges(th, dom, r1[0], 1e3, g.nbr, g.dist, g.axis, g.area), (1, 4))
    assert np.array_equal(d1, d4)


def test_chart_table_matches_at_rounding_boundaries_and_face_ties():
    for n in (48, 217):
        half_cells = np.arange(n)+.5
        ratios = np.tan(cube.QUARTER*(2*half_cells/n-1))
        trials = np.concatenate((ratios-1e-12, np.nextafter(ratios, -np.inf), ratios,
                                 np.nextafter(ratios, np.inf), ratios+1e-12, [-1., 0., 1.]))[:, None]
        # Vary each chart axis at the half-cell thresholds and at face edges, including three-way ties.
        p = np.concatenate([normal + trials*axis + edge*other
                            for normal, right, down in zip(cube.NORMAL, cube.RIGHT, cube.DOWN)
                            for axis, other in ((right, down), (down, right))
                            for edge in (-1., 0., 1.)] + [np.zeros((1, 3))])
        p = np.concatenate((p, p*1e-300, p*1e300))
        args = (p, n, cube.NORMAL, cube.RIGHT, cube.DOWN, fields._FACE_POS, fields._FACE_NEG)
        assert np.array_equal(fields.nearest_entries(*args, fields._chart_table(n)), fields.nearest_entries(*args))


def test_forward_map_and_disc_caches_return_what_they_replace():
    g = CubeSphere(32)
    sim = PlateSim(g, seed=5, t0=-4.)
    for _ in range(2):
        sim.step()
    for pl in sim.plates:
        wn = sim._forward_nodes(pl)
        assert np.array_equal(wn, fields.nearest_nodes(g, g.xyz@pl.Q.T)) and not wn.flags.writeable
    m = sim.mantle
    h = m.hotspots[0]
    w = m._disc(h['pos'], 5e5)
    ang = np.arccos(np.clip(g.xyz@h['pos'], -1, 1))*g.radius
    assert np.array_equal(w, np.exp(-.5*(ang/5e5)**2)) and m._disc(h['pos'], 5e5) is w
    m.hotspots = m.hotspots[1:]
    m.advance(sim.t, 0.)                               # a dead hotspot's discs leave the cache
    assert all(k[0] != h['pos'].tobytes() for k in m._discs)
