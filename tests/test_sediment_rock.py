"""sediment_rock (pass 95, journal 126): the sediment layer's lithified part. 'none' (as built): all of it loose;
'burial': sediment buried below rock_depth_m (bulk) is rock, the start's continental cover is rock, and only the loose
top is the cover the 'cover' law entrains without bedrock's threshold."""
import numpy as np
import pytest

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.isostasy import bulk_thickness, loose_cap, solid_above


def test_default_is_burial_and_none_keeps_the_rock_layer_empty():
    assert P.PARAMS['sediment_rock'] == 'burial'           # the Rock Cover Sweep (journal 129)
    sim = P.PlateSim(CubeSphere(48), seed=3, params={'threads': 1, 'sediment_rock': 'none'})
    sim.state()
    for _ in range(4):
        sim.step()
        sim.state()
    assert not sim.world['rock'].any()
    assert all(not pl.rock.any() for pl in sim.plates)


def test_unknown_rule_raises_and_burial_needs_the_cover_law():
    with pytest.raises(ValueError):
        P.PlateSim(CubeSphere(16), seed=1, params={'sediment_rock': 'nope', 'threads': 1})
    sim = P.PlateSim(CubeSphere(16), seed=1, params={'sediment_rock': 'burial', 'sediment_entrainment': 'bedrock',
                                                      'threads': 1})
    sim.state()
    with pytest.raises(ValueError):
        sim.step()


def test_loose_cap_is_the_solid_in_the_top_rock_depth():
    d = np.array([10., 500., 2000., 6000.])
    for phi0, c in ((.56, .39), (.63, .51)):
        s = solid_above(d, phi0, c)
        assert np.allclose(bulk_thickness(s, phi0, c), d, atol=1e-3)
    cap = loose_cap(np.array([True, False]), P.PARAMS)
    assert 1150 < cap[0] < 1300 and 1150 < cap[1] < 1300   # ~1.22 km of solid in 2 km of bulk


def test_burial_start_cover_is_rock_and_the_rest_of_the_start_unchanged():
    g = CubeSphere(48)
    a = P.PlateSim(g, seed=3, params={'threads': 1})
    b = P.PlateSim(g, seed=3, params={'threads': 1, 'sediment_rock': 'burial'})
    for k in ('cont', 'thick', 'born', 'orog', 'therm', 'sed'):
        assert np.array_equal(a.world[k], b.world[k], equal_nan=True), k
    assert np.array_equal(a.owner, b.owner)
    w = b.world; c = w['cont']
    assert w['sed'][c].sum() > 0
    assert np.allclose(w['rock'][c], w['sed'][c], rtol=1e-6)            # Earth's stable cover is rock
    cap = loose_cap(c, b.p)
    assert np.allclose(w['rock'][~c], np.maximum(w['sed'][~c]-cap[~c], 0.), atol=1e-2)


def test_burial_keeps_rock_within_the_sediment_and_lithifies_below_the_cap():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=5, params={'threads': 1, 'sediment_rock': 'burial'})
    sim.state()
    for _ in range(12):
        sim.step()
        w = sim.world
        sed = np.maximum(w['sed'], 0.)
        assert (w['rock'] >= 0).all() and (w['rock'] <= sed+1e-6).all()
        assert (w['rock'] >= sed-loose_cap(w['cont'], sim.p)-1e-3).all()
        for pl in sim.plates:
            m = pl.mask
            assert (pl.rock[m] >= 0).all() and (pl.rock[m] <= np.maximum(pl.sed[m], 0.)+1e-3).all()
        sim.state()


def test_burial_keeps_the_ledgers_and_holds_the_cover_on_the_land():
    g = CubeSphere(48)
    out = {}
    for rule in ('none', 'burial'):
        sim = P.PlateSim(g, seed=5, params={'threads': 1, 'sediment_rock': rule})
        sim.state()
        v0, s0 = sim.continental_volume(), sim.continental_sediment_volume()
        for _ in range(12):
            sim.step()
            sim.state()
        assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
        s1 = sim.continental_sediment_volume()
        assert abs((s1-s0)-sum(sim.sed_ledger.values())) < 1e-6*max(s1, 1.)+1e-3
        out[rule] = sim.erosion_totals
    # rock does not leave below bedrock's threshold: the land loses less of its surface
    assert -out['burial']['continental_eroded_m3'] < -out['none']['continental_eroded_m3']
