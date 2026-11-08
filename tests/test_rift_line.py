"""Pass 71 (journal 94): `rift_line`, where a split's line runs. 'voronoi' (as built) and 'weak' (the line lies in weak
lithosphere: C's thermal thickness and the belts; Audet & Buergmann 2011; the default since pass 72). Pass 72 (journal
95): 'chord', the least-force cut through a continental nucleus that leaves a plate on both sides (the default since
pass 73). Pass 73 (journal 96): 'force', the candidate cut closest to failure (opening force over strength)."""
import numpy as np
import pytest

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P


def _sim(params, seed=5):
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=seed, params={'water_volume_km3': None, **params})
    sim.step()
    return sim


def _line_nodes(sim, own, lab):
    nb = sim.g.nbr; ok = (nb >= 0) & own[np.maximum(nb, 0)]
    return own & np.any(ok & (lab[np.maximum(nb, 0)] != lab[:, None]), 1)


def _draws(sim, n=40):
    pl = max(sim.plates, key=lambda q: (sim.owner == q.id).sum())
    own = sim.owner == pl.id
    sim.rng = np.random.default_rng(11)
    return own, [_line_nodes(sim, own, sim._rift_line(pl, own)) for _ in range(n)]


def test_default_is_force_and_unknown_raises():
    assert P.PARAMS['rift_line'] == 'force'
    with pytest.raises(ValueError):
        _sim({'rift_line': 'nope'})


def test_weak_line_lies_in_young_lithosphere():
    """A band of young (thin) continental lithosphere across the plate: 'weak' lines run in it more than 'voronoi'."""
    shares = {}
    for rule in ('voronoi', 'weak'):
        sim = _sim({'rift_line': rule, 'lithosphere': 'thermal_age'})
        w = sim.world
        w['cont'][:] = True
        band = np.abs(sim.g.xyz[:, 0]) < .15
        w['therm'][:] = sim.t-2500.                    # old everywhere (L 194 km) ...
        w['therm'][band] = sim.t                       # ... but a young band (L 94 km)
        w['orog'][:] = np.nan
        own, lines = _draws(sim)
        shares[rule] = np.array([(ln & band).sum()/max(ln.sum(), 1) for ln in lines])
    _assert_more(shares)


def _assert_more(shares):
    """Paired draws (same seeds and wiggle; only the cost differs): 'weak' puts more of its line in the band."""
    v, k = shares['voronoi'], shares['weak']
    assert k.mean() > 1.3*v.mean() and (k >= v).mean() >= .75, (v.mean(), k.mean(), (k >= v).mean())


def test_weak_line_follows_belts_voronoi_avoids_them():
    """The belt factor: 'weak' divides by it (belts slow both fronts), 'voronoi' multiplies (belts let one through)."""
    shares = {}
    for rule in ('voronoi', 'weak'):
        sim = _sim({'rift_line': rule, 'lithosphere': 'none'})
        w = sim.world
        band = np.abs(sim.g.xyz[:, 1]) < .15
        w['orog'][:] = np.nan
        w['orog'][band] = sim.t
        own, lines = _draws(sim)
        shares[rule] = np.array([(ln & band).sum()/max(ln.sum(), 1) for ln in lines])
    _assert_more(shares)


def test_strength_is_c_thickness_over_earth_mean():
    sim = _sim({'rift_line': 'weak', 'lithosphere': 'thermal_age'})
    w = sim.world
    st = sim._line_strength()
    c = w['cont'] & np.isfinite(sim.t-w['therm'])
    p = sim.p
    exp = (p['lith_L0_km']+p['lith_rate_km_myr']*np.maximum(sim.t-w['therm'][c], 0.))/p['lith_ref_km']
    assert np.allclose(st[c], exp) and np.all(st[~c] == 1.)
    assert np.all(_sim({'rift_line': 'weak', 'lithosphere': 'none'})._line_strength() == 1.)


def test_weak_line_splits_and_keeps_the_ledger():
    from test_rift_speed import _Always
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'rift_line': 'weak', 'split_min_area_mkm2': 1.})
    v0 = sim.continental_volume()
    sim.step()
    sim.t_start = sim.t-90.                              # Earth's breakup peak (160 Ma)
    rng = sim.rng
    sim.rng = _Always(rng)
    sim._split_events_earth()                            # every candidate split taken
    sim.rng = rng
    assert any(e['kind'] == 'rift' for e in sim.events)
    for _ in range(20):
        sim.step()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0


def _big_continental_plate(rule, seed=5):
    sim = _sim({'rift_line': rule, 'lithosphere': 'thermal_age'}, seed=seed)
    pl = max(sim.plates, key=lambda q: (sim.owner == q.id).sum())
    return sim, pl, sim.owner == pl.id


def test_chord_leaves_a_plate_on_both_sides_and_crosses_the_continent():
    sim, pl, own = _big_continental_plate('chord')
    g = sim.g
    floor = sim.p['piece_min_area_mkm2']*1e12
    sim.rng = np.random.default_rng(11)
    n0 = sum(getattr(sim, 'chord_stats', {}).values())
    kinds = []
    for _ in range(15):
        lab = sim._rift_line(pl, own)
        kinds.append(sim._line_kind)
        if sim._line_kind != 'chord':
            continue
        part = own & (lab == 1)
        assert g.area[part].sum() >= floor and g.area[own & ~part].sum() >= floor
        assert sim._split_margin_km(own, lab) > 0
        from planet_tectonics.graphops import components
        assert components(part, g.nbr).max() == 0       # the new part is one piece
    assert kinds.count('chord') >= 12, kinds
    assert sum(sim.chord_stats.values())-n0 == 15


def test_chord_parts_are_smaller_than_the_bisector_halves():
    """The bisector halves its plate; the chord cuts where the least force does, often a smaller piece."""
    areas = {}
    for rule in ('weak', 'chord'):
        sim, pl, own = _big_continental_plate(rule)
        sim.rng = np.random.default_rng(11)
        a = []
        for _ in range(20):
            part = own & (sim._rift_line(pl, own) == 1)
            a.append(min(sim.g.area[part].sum(), sim.g.area[own & ~part].sum()))
        areas[rule] = np.median(a)
    assert areas['chord'] < .5*areas['weak'], areas


def test_chord_falls_back_to_the_bisector_and_ocean_splits_keep_it():
    sim, pl, own = _big_continental_plate('chord')
    sim.world['cont'][:] = False
    lab = sim._rift_line(pl, own)
    assert sim._line_kind == 'bisector' and set(np.unique(lab[own])) == {0, 1}
    n = dict(sim.chord_stats)
    sim._rift_line(pl, own, chord=False)                 # the ocean-only branch draws the bisector, uncounted
    assert sim.chord_stats == n


def test_chord_split_keeps_the_ledger():
    from test_rift_speed import _Always
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'rift_line': 'chord', 'split_min_area_mkm2': 1.})
    v0 = sim.continental_volume()
    sim.step()
    sim.t_start = sim.t-90.
    rng = sim.rng
    sim.rng = _Always(rng)
    sim._split_events_earth()
    sim.rng = rng
    rifts = [e for e in sim.events if e['kind'] == 'rift' and 'line' in e]
    assert rifts and all(e['line'] in ('chord', 'bisector') for e in rifts)
    for _ in range(20):
        sim.step()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0


def _forced_continental_plate():
    """A plate with continent and slab pull: the test supercontinent overrides all round (no boundary force under the
    fit's zero suction), so its rim is handed a trench of its own (the downgoing side) on a third of its edge."""
    sim, pl, own = _big_continental_plate('force')
    g = sim.g
    edge = own & np.any((g.nbr >= 0) & ~own[np.clip(g.nbr, 0, None)], 1)
    c = g.xyz[own].mean(0); c /= np.linalg.norm(c)
    k = np.cross(c, [0., 0., 1.]); k /= np.linalg.norm(k)
    side = edge & (g.xyz@k > np.quantile(g.xyz[edge]@k, 2/3))
    sim.trench = sim.trench | side
    sim.btype = np.where(side, 2, sim.btype).astype(sim.btype.dtype)
    sim._cut_forces_t = None
    return sim, pl, own


def test_force_cut_leaves_plates_and_is_in_tension():
    sim, pl, own = _forced_continental_plate()
    g = sim.g
    floor = sim.p['piece_min_area_mkm2']*1e12
    sim.rng = np.random.default_rng(11)
    kinds = []
    from planet_tectonics.graphops import components
    for _ in range(15):
        lab = sim._rift_line(pl, own)
        kinds.append(sim._line_kind)
        if sim._line_kind != 'force':
            continue
        part = own & (lab == 1)
        assert g.area[part].sum() >= floor and g.area[own & ~part].sum() >= floor
        assert components(part, g.nbr).max() == 0
        assert sim._cut_tension(own, part) > 0           # in tension, read the same from either side
        assert np.isclose(sim._cut_tension(own, part), sim._cut_tension(own, own & ~part), rtol=1e-6)
    assert kinds.count('force') >= 10, kinds
    assert sim.chord_stats.get('force', 0) >= 10


def test_force_cut_is_nearer_failure_than_the_cheapest():
    """Same nuclei and candidates: the force cut's opening force over its strength is not below the cheapest cut's."""
    r = {}
    for rule in ('chord', 'force'):
        sim, pl, own = _forced_continental_plate()
        sim.p['rift_line'] = rule
        sim.rng = np.random.default_rng(11)
        out = []
        for _ in range(12):
            lab = sim._rift_line(pl, own)
            if sim._line_kind == 'bisector':
                continue
            part = own & (lab == 1)
            out.append(sim._cut_tension(own, part)/max(sim._split_margin_km(own, lab), 1.))
        r[rule] = np.median(out)
    assert r['force'] > r['chord'], r


def test_force_without_tension_takes_the_cheapest_cut():
    """No boundary force (the test supercontinent overrides all round): the cut is the chord's, same rng."""
    labs = {}
    for rule in ('chord', 'force'):
        sim, pl, own = _big_continental_plate(rule)
        sim.rng = np.random.default_rng(11)
        labs[rule] = [sim._rift_line(pl, own) for _ in range(5)]
        if rule == 'force':
            assert set(sim.chord_stats) <= {'chord', 'bisector'}
    assert all(np.array_equal(a, b) for a, b in zip(labs['chord'], labs['force']))


def test_force_split_keeps_the_ledger():
    from test_rift_speed import _Always
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'rift_line': 'force', 'split_min_area_mkm2': 1.})
    v0 = sim.continental_volume()
    sim.step()
    sim.t_start = sim.t-90.
    rng = sim.rng
    sim.rng = _Always(rng)
    sim._split_events_earth()
    sim.rng = rng
    rifts = [e for e in sim.events if e['kind'] == 'rift' and 'line' in e]
    assert rifts and all(e['line'] in ('force', 'chord', 'bisector') for e in rifts)
    for _ in range(20):
        sim.step()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
