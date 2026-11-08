"""Pass 71 (journal 94): `rift_line`, where a split's line runs. 'voronoi' (as built, the default) and 'weak' (the line
lies in weak lithosphere: C's thermal thickness and the belts; Audet & Buergmann 2011)."""
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


def test_default_is_voronoi_and_unknown_raises():
    assert P.PARAMS['rift_line'] == 'voronoi'
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
