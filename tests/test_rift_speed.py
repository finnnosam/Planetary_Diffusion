"""Pass 70 (journal 93): `rift_speed_cm_yr`, the opening of a split that crosses a continent (Earth's fast rift
phase, Brune et al. 2016). Off (0) by default: such splits open at `split_speed_cm_yr` as built."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P


class _Always:
    """The run's rng with uniform() always 0, so every candidate split is taken; the line is still drawn as usual."""
    def __init__(self, rng):
        self._rng = rng

    def uniform(self, *a, **k):
        return 0.

    def __getattr__(self, name):
        return getattr(self._rng, name)


def _speeds(params):
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=5, params={'water_volume_km3': None, **params})
    sim.step()
    sim.t_start = sim.t-90.                              # Earth's breakup peak (160 Ma)
    got = []
    real = sim.split
    sim.split = lambda pl, part, speed, kind: got.append((kind, speed, sim._split_margin_km(sim.owner == pl.id,
                                                          part.astype(np.int64)))) or real(pl, part, speed, kind)
    sim.rng = _Always(sim.rng)
    sim._split_events_earth()
    return sim, got


def test_rift_speed_is_off_by_default():
    assert P.PARAMS['rift_speed_cm_yr'] == 0.
    sim, got = _speeds({})
    assert got and all(s == P.PARAMS['split_speed_cm_yr'] for _, s, _ in got)
    assert sim._rift_speed(500.) == sim._rift_speed(0.) == P.PARAMS['split_speed_cm_yr']


def test_a_continental_split_opens_at_the_rift_speed():
    sim, got = _speeds({'rift_speed_cm_yr': 3.5})
    cont = [s for _, s, L in got if L > 0]
    assert cont and all(s == 3.5 for s in cont)
    assert sim._rift_speed(0.) == P.PARAMS['split_speed_cm_yr']      # ocean-only splits keep split_speed
    # the halves' Euler vectors differ by exactly the rift speed at the split
    pl = max(sim.plates, key=lambda q: (sim.owner == q.id).sum())
    own = sim.owner == pl.id
    half = own & (sim.g.xyz@sim.g.xyz[np.flatnonzero(own)[0]] > np.median(sim.g.xyz[own]@sim.g.xyz[np.flatnonzero(own)[0]]))
    w0 = sim._omega_vec(pl)
    new = P.PlateSim.split(sim, pl, half, 3.5, 'test')
    assert new is not None and abs(np.linalg.norm(sim._omega_vec(new)-sim._omega_vec(pl))-3.5) < 1e-9
    assert np.allclose(.5*(sim._omega_vec(new)+sim._omega_vec(pl)), w0) or min(new.speed, pl.speed) == 1e-3


def test_rift_speed_keeps_the_ledger():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'rift_speed_cm_yr': 3.5, 'split_min_area_mkm2': 1.})
    v0 = sim.continental_volume()
    for _ in range(30):
        sim.step()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
