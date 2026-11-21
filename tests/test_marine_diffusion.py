"""Pass 97 (journal 130): marine_diffusion. 'uniform' (as built) runs the landscape kernel's marine diffusion
(Kd_marine, at every depth to 1,000 m) on the plate grid; 'none' turns it off there, leaving the shelf law (wave
transport decaying with depth) and the gravity paths (deltas, graded shelves, the margin wedge, turbidites)."""
import numpy as np
import pytest

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P


def _sim(rule, **kw):
    sim = P.PlateSim(CubeSphere(48), seed=3, params={'threads': 1, 'marine_diffusion': rule, **kw})
    sim.state()
    return sim


def test_default_is_uniform_and_unknown_rules_raise():
    assert P.PARAMS['marine_diffusion'] == 'uniform'
    sim = _sim('wave')
    with pytest.raises(ValueError):
        sim._erode(sim.world['thick'].astype(float))


def test_none_changes_the_erosion():
    out = {}
    for rule in ('uniform', 'none'):
        sim = _sim(rule)
        th = sim.world['thick'].astype(float)
        e, ds = sim._erode(th)
        out[rule] = (e, ds)
    assert not np.array_equal(out['uniform'][1], out['none'][1])


def test_none_does_not_cut_deep_submerged_continent():
    """Below the shelf law's reach and off the gravity paths, nothing removes material without the uniform term: no
    node deeper than 700 m (with every neighbour deeper than 650 m) loses crust or sediment under 'none'."""
    cut = {}
    for rule in ('uniform', 'none'):
        sim = _sim(rule)
        th = sim.world['thick'].astype(float)
        z0 = sim.elevation(th, sim.world['sed'])
        nb = sim.g.nbr
        deep = (z0 < -700.) & np.all(np.where(nb >= 0, z0[np.clip(nb, 0, None)] < -650., True), axis=1)
        e, ds = sim._erode(th)
        cut[rule] = float(np.minimum(e+ds, 0.)[deep].sum())
    assert cut['none'] == 0.
    assert cut['uniform'] < 0.


def test_none_keeps_the_ledgers():
    sim = _sim('none')
    v0, s0 = sim.continental_volume(), sim.continental_sediment_volume()
    for _ in range(10):
        sim.step()
        sim.state()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
    s1 = sim.continental_sediment_volume()
    assert abs((s1-s0)-sum(sim.sed_ledger.values())) < 1e-6*max(s1, 1.)+1e-3
    st = sim.state()
    assert np.isfinite(st['elevation_m']).all() and (st['sediment_m'] >= 0).all()
