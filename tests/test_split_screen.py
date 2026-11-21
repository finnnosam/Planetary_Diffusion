"""The split candidates' screen (journal 133): graphops.two_big_pieces must never reject a cut that the exact labelling
accepts, and _chord_line must give the same cut with and without it."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import events as ev
from planet_tectonics import plates as P
from planet_tectonics.graphops import components, shortest_paths, trace_cut, two_big_pieces


def _exact_two_big(blocked, lnbr, area, floor):
    comp = components(~blocked, lnbr)
    n = comp.max()+1
    if n < 2:
        return False
    a = np.sort(np.bincount(comp[comp >= 0], weights=area[comp >= 0], minlength=n))
    return bool(a[-2] >= floor)


def test_the_screen_never_rejects_a_cut_the_exact_count_accepts():
    g = CubeSphere(24)
    rng = np.random.default_rng(5)
    seen = {True: 0, False: 0}
    for trial in range(60):
        c = rng.normal(size=3); c /= np.linalg.norm(c)
        own = g.xyz@c > np.cos(rng.uniform(.4, 1.2))
        if trial % 3 == 0:                               # a plate with detached crumbs, as real plates carry
            own |= rng.uniform(size=g.size) < .01
        idx = np.flatnonzero(own)
        loc = np.full(g.size, -1, np.int64); loc[idx] = np.arange(idx.size)
        nb = g.nbr[idx]
        lnbr = np.where(nb >= 0, loc[np.clip(nb, 0, None)], -1)
        cost = np.where((g.nbr >= 0) & own[:, None] & own[np.clip(g.nbr, 0, None)], g.dist*rng.uniform(.5, 2, g.nbr.shape), 1e30)
        area = g.area[idx]
        comp0 = components(np.ones(idx.size, bool), lnbr)
        for _ in range(10):
            s = int(rng.choice(idx))
            d, pred = shortest_paths(s, g.nbr, cost)
            reach = idx[d[idx] < 1e29]
            cut = np.zeros(g.size, bool)
            for e in rng.choice(reach, 2):
                trace_cut(cut, pred, s, int(e), g.nbr, g.axis)
            bl = cut[idx]
            floor = rng.uniform(.02, .4)*area.sum()
            big0 = np.flatnonzero(np.bincount(comp0, weights=area) >= floor*(1-1e-9))
            if not (big0.size == 0 or (big0.size == 1 and comp0[loc[s]] == big0[0])):
                continue                                 # _chord_line does not screen these
            exact = _exact_two_big(bl, lnbr, area, floor)
            screen = two_big_pieces(bl, lnbr, area, floor)
            assert screen or not exact
            seen[exact] += 1
    assert seen[True] > 5 and seen[False] > 5


def test_chord_line_is_the_same_with_and_without_the_screen(monkeypatch):
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=11, params={'water_volume_km3': None, 'threads': 1})
    for _ in range(4):
        sim.step(); sim.state()
    plates = [pl for pl in sim.plates if (sim.owner == pl.id).any()]
    assert plates
    for pl in plates:
        own = sim.owner == pl.id
        state = sim.rng.bit_generator.state
        a = sim._chord_line(own)
        sim.rng.bit_generator.state = state
        monkeypatch.setattr(ev, 'two_big_pieces', lambda *args: True)
        b = sim._chord_line(own)
        monkeypatch.undo()
        assert (a is None and b is None) or np.array_equal(a, b)
