"""backarc_extension (ninety-second pass, journal 119): the extensional branch of the upper-plate rate u - rollback.
Volume moves from the back-arc band to the node the arc-forearc block steps onto; the ledgers close; on by default
since pass 93 (journal 122), off books nothing."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P

BASE = {'water_volume_km3': None, 'start_sediment': 'none'}


def test_backarc_default_and_off_books_nothing():
    assert P.PARAMS['backarc_extension'] is True and P.PARAMS['rollback_mm_yr'] == 9.
    sim = P.PlateSim(CubeSphere(48), seed=3, params={**BASE, 'backarc_extension': False})
    for _ in range(6):
        sim.step()
    assert 'backarc' not in sim.ledger and sim.erosion_totals.get('backarc_steps', 0) == 0
    assert sim.erosion_totals['overriding_trench_m'] > 0       # the regime record runs in every run (output only)


def test_backarc_run_moves_volume_and_keeps_the_ledgers():
    # rollback 200 mm/yr: every continental overriding segment extends, so the mechanism acts on the coarse grid
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={**BASE, 'backarc_extension': True, 'rollback_mm_yr': 200.})
    v0 = sim.continental_volume()
    s0 = sim.continental_sediment_volume()
    thin_seen = False
    for _ in range(20):
        sim.step()
        thin_seen |= bool((sim.tectonic < 0).any())
    ET = sim.erosion_totals
    assert ET['backarc_steps'] > 0 and ET['backarc_stretched_m3'] > 0 and thin_seen
    assert sim.ledger['backarc'] > 0 and sim.area_ledger['backarc'] > 0
    # the strips carry what the bands gave (each step's volume is booked once, on the plate raster)
    assert abs(sim.ledger['backarc']-ET['backarc_stretched_m3']) < 2e-2*ET['backarc_stretched_m3']
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
    s1 = sim.continental_sediment_volume()
    assert abs((s1-s0)-sum(sim.sed_ledger.values())) < 1e-6*max(s1, 1.)+1e-3


def test_backarc_band_follows_earths_profile():
    # the band weight is Earth's: nothing within 200 km of the trench or beyond 1,500 km
    from planet_tectonics.params import BACKARC_PROFILE_KM, BACKARC_PROFILE_W
    km, w = np.asarray(BACKARC_PROFILE_KM), np.asarray(BACKARC_PROFILE_W)
    assert km[0] == 200. and km[-1] == 1500. and len(km) == len(w)+1 and np.all(w > 0)
    mids = .5*(km[1:]+km[:-1])
    mu = (w*mids).sum()/w.sum()
    assert 700 < mu < 790                                     # excess-weighted mean 745 km (thin_by_trench.txt)


def test_backarc_needs_margin_shortening():
    import pytest
    sim = P.PlateSim(CubeSphere(48), seed=3, params={**BASE, 'backarc_extension': True, 'margin_shortening': False})
    with pytest.raises(ValueError):
        for _ in range(3):
            sim.step()


def test_lallemand_rollback_follows_the_subducting_plate():
    # rollback_rule 'lallemand': the trench rolls back at 14 - 0.24 Vsub mm/yr (SB04); the rule changes the rate only
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={**BASE, 'rollback_rule': 'lallemand'})
    for _ in range(4):
        sim.step()
    src = (sim.btype == P.SUBDUCTION) & ~sim.trench
    idx = np.flatnonzero(src)
    nb = g.nbr[idx]
    other = (nb >= 0) & (sim.owner[np.maximum(nb, 0)] != sim.owner[idx][:, None])
    u, n = sim._trenchward_speed(idx, nb, other)
    r = sim._rollback(idx, nb, other, n)
    assert r.shape == idx.shape and np.isfinite(r).all() and np.ptp(r) > 0
    # one node by hand: the mean landward speed of its oceanic neighbours on other plates
    i = 0
    x = g.xyz[idx[i]]
    js = [j for j, o in zip(nb[i], other[i]) if o and not sim.world['cont'][j]] or [j for j, o in zip(nb[i], other[i]) if o]
    om = {pl.id: pl.omega(g.radius) for pl in sim.plates}
    vs = np.mean([-np.dot(np.cross(om[int(sim.owner[j])], x)*g.radius, n[i]) for j in js])
    assert abs(r[i]-(14e-3-.24*vs)) < 1e-12
    sim.p['rollback_rule'] = 'fixed'
    assert sim._rollback(idx, nb, other, n) == sim.p['rollback_mm_yr']*1e-3
