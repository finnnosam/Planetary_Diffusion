"""Pass 77 (journal 101): collision_lithosphere 'thicken', the lithosphere thickening with collisional crust."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.lithosphere import thickened_therm, thickness_km


def test_default_is_reset():
    assert P.PARAMS['collision_lithosphere'] == 'reset'


def test_pure_shear_thickens_L_by_the_crust_factor():
    p = P.PARAMS
    t = -10.
    therm = np.array([t-1000., t-0., t-500., np.nan, t-2000.])
    thick = np.array([40e3, 40e3, 35e3, 40e3, 40e3])
    gain = np.array([10e3, 8e3, 0., 5e3, -1e3])
    new, m = thickened_therm(therm, thick, gain, t, p)
    assert list(m) == [True, True, False, False, False]
    for i in (0, 1):
        f = (thick[i]+gain[i])/thick[i]
        assert abs(thickness_km(t-new[i], p)-thickness_km(t-therm[i], p)*f) < 1e-9
    assert new[2] == therm[2] and np.isnan(new[3]) and new[4] == therm[4]


def test_thicken_resets_only_the_arc_belt():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'collision_lithosphere': 'thicken', 'threads': 1})
    sim.state()
    cont = sim.world['cont']
    collided = cont.copy()
    sim.belt = np.zeros(g.size, bool); sim.belt[np.flatnonzero(cont)[:5]] = True
    r = sim._thermal_resets(collided)
    assert r.sum() == 5 and np.all(r == (sim.belt & cont))
    sim.p = {**sim.p, 'collision_lithosphere': 'reset'}
    assert sim._thermal_resets(collided).sum() >= cont.sum()


def test_thicken_runs_and_keeps_thick_crust_on_thicker_lithosphere():
    g = CubeSphere(48)
    out = {}
    for rule in ('reset', 'thicken'):
        sim = P.PlateSim(g, seed=5, params={'collision_lithosphere': rule, 'threads': 1})
        sim.state()
        for _ in range(40):
            sim.step()
            st = sim.state()
        out[rule] = st
    def L_thick(st):
        T = st['crust_thickness_m']; m = st['continental'] & (T > 42e3) & np.isfinite(st['lithosphere_age_myr'])
        return thickness_km(st['lithosphere_age_myr'][m], P.PARAMS).mean() if m.any() else np.nan
    assert L_thick(out['thicken']) > L_thick(out['reset'])
