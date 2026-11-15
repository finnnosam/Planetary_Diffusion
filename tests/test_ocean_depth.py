"""Pass 82 (journal 107): ocean_depth_rule, the ocean floor's basement depth by age. 'gdh1' (Stein & Stein 1992; the
default to pass 82) or Earth's floor measured on the n = 217 graph ('earth': normal oceanic crust; 'earth_all': every dated floor)."""
import numpy as np
import pytest

from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import isostasy as iso
from planet_tectonics import plates as P

REF = OUTPUTS/'planet-reference'
HAVE = all((REF/f).exists() for f in ('earth-n217.npz', 'earth-seafloor-age-n217.npz', 'globsed-n217.npz')) and (
    DATA/'crust1'/'CNtype1-1.txt').exists()


def test_default_is_earth():
    """'earth' is the default since pass 83 (journal 108); 'gdh1' still gives GDH1 exactly."""
    assert P.PARAMS['ocean_depth_rule'] == 'earth'
    a = np.array([0., 1., 19.9, 20., 80., 250.])
    assert np.array_equal(P.ocean_depth(a), P.ocean_depth(a, {**P.PARAMS, 'ocean_depth_rule': 'earth'}))
    assert np.array_equal(P.ocean_depth(a, {**P.PARAMS, 'ocean_depth_rule': 'gdh1'}), P.gdh1_depth(a))


@pytest.mark.skipif(not HAVE, reason='Earth references on n = 217 (ETOPO, Seton ages, GlobSed) or CRUST1.0 not present')
@pytest.mark.parametrize('rule', ['earth', 'earth_all'])
def test_table_reproduces_earth(rule):
    t, d = iso.earth_floor_table(rule)
    assert np.allclose(t, iso.EARTH_FLOOR_AGE_MYR[rule], atol=6e-3)
    assert np.allclose(d, iso.EARTH_FLOOR_DEPTH_M[rule], atol=.06)


@pytest.mark.parametrize('rule', ['earth', 'earth_all'])
def test_law_is_monotone_continuous_and_passes_through_its_knots(rule):
    p = {**P.PARAMS, 'ocean_depth_rule': rule}
    T = np.array(iso.EARTH_FLOOR_AGE_MYR[rule]); D = np.array(iso.EARTH_FLOOR_DEPTH_M[rule])
    assert len(T) == len(D) and np.all(np.diff(T) > 0) and np.all(np.diff(D) >= 0)
    assert np.allclose(P.ocean_depth(T, p), D)
    a = np.linspace(0., 320., 32001)
    d = P.ocean_depth(a, p)
    assert np.all(np.diff(d) >= -1e-9)                                  # the floor never rises as it ages
    assert np.abs(np.diff(d))[a[1:] > .5].max() < 3.                    # no step (0.01 Myr apart) once past sqrt(age)'s start
    assert abs(P.ocean_depth(np.array([T[0]+1e-6]), p)[0]-P.ocean_depth(np.array([T[0]-1e-6]), p)[0]) < 1e-2   # nor at the first knot
    assert d[-1] == D[-1] and P.ocean_depth(np.array([-5.]), p)[0] == d[0]   # flat beyond the last knot; age < 0 as 0
    # below the first knot: GDH1's own shape, shifted through the knot
    assert abs((D[0]-d[0])-(P.gdh1_depth(T[0])-P.gdh1_depth(0.))) < 1e-9


def test_earths_floor_stands_above_gdh1_beyond_10_myr():
    a = np.array([15., 30., 60., 100., 150., 250.])
    g = P.gdh1_depth(a)
    e = P.ocean_depth(a, {**P.PARAMS, 'ocean_depth_rule': 'earth'})
    al = P.ocean_depth(a, {**P.PARAMS, 'ocean_depth_rule': 'earth_all'})
    assert np.all((g-e > 100.) & (g-e < 260.))                          # normal floor: 120-240 m above GDH1
    assert np.all(al <= e+1e-9) and np.all(g-al < 400.)                 # with the plateaus in the mean, shallower still
    assert abs(P.ocean_depth(np.array([300.]), {**P.PARAMS, 'ocean_depth_rule': 'earth'})[0]-5413.) < .1


def test_unknown_rule_raises():
    with pytest.raises(ValueError):
        P.ocean_depth(np.array([10.]), {**P.PARAMS, 'ocean_depth_rule': 'nope'})


def test_pooling_keeps_the_weighted_mean_and_orders_the_sequence():
    v = np.array([1., 3., 2., 2.5, 5., 4.]); w = np.array([1., 2., 1., 1., 3., 1.])
    q = iso._pool_increasing(v, w)
    assert np.all(np.diff(q) >= 0) and abs((q*w).sum()-(v*w).sum()) < 1e-12
    assert np.allclose(q, [1., 2.625, 2.625, 2.625, 4.75, 4.75])
    assert np.array_equal(iso._pool_increasing([1., 2., 3.], [1., 1., 1.]), [1., 2., 3.])


def test_isostatic_elevation_and_the_surface_follow_the_rule():
    cont = np.array([False, False]); thick = np.array([7e3, 7e3]); age = np.array([0., 100.])
    for rule in ('gdh1', 'earth'):
        p = {**P.PARAMS, 'ocean_depth_rule': rule}
        assert np.allclose(P.isostatic_elevation(cont, thick, age, p), -P.ocean_depth(age, p))
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'threads': 1})
    for _ in range(3):
        sim.step()
    w = sim.world
    age = np.where(w['cont'], 0., sim.t-w['born'])
    plain = (~w['cont']) & (np.abs(w['thick']-sim.p['ocean_thick_m']) < 1.) & (w['sed'] < 50.) & (age > 30.)
    sh = P.sediment_column(w['sed'], w['cont'], sim.p)[0]
    assert plain.sum() > 200
    out = {}
    for rule in ('gdh1', 'earth'):
        sim.p = {**sim.p, 'ocean_depth_rule': rule}
        z = sim.elevation()                                             # relative to the sea it solves
        out[rule] = (z+sim.sea_level, sim.sea_level)
        own = -P.ocean_depth(age, sim.p)+sh+sim.dyn                     # a plain floor node's own level on the datum
        assert abs(np.median((out[rule][0]-own)[plain])) < 30.          # its thin sediment's load and far flexure only
    lift = (out['earth'][0]-out['gdh1'][0])[plain]
    want = (P.gdh1_depth(age)-P.ocean_depth(age, {**sim.p, 'ocean_depth_rule': 'earth'}))[plain]
    assert np.allclose(lift, want, atol=5.)                             # the floor rises by the laws' difference
    assert out['earth'][1] > out['gdh1'][1]+50.                         # and the sea with it


def test_a_run_on_earths_floor_keeps_the_young_ocean_gate():
    from planet_reviews import elevation_check as ec
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=5, params={'threads': 1, 'ocean_depth_rule': 'earth'})
    sim.state()
    for _ in range(25):
        sim.step(); st = sim.state()
    z = st['elevation_m']
    assert np.isfinite(z).all() and z.max() < 9000.
    w = sim.world
    age = np.where(w['cont'], np.nan, sim.t-w['born'])
    young = (~w['cont']) & (age < 20) & (np.abs(w['thick']-sim.p['ocean_thick_m']) < 1000.) & (w['sed'] < 300.)
    own = -P.ocean_depth(np.nan_to_num(age), sim.p)+(w['thick']-sim.p['ocean_thick_m'])+sim.dyn-sim.sea_level
    assert abs(np.median((z-own)[young])) < 300.                        # the gate's young-ocean offset, on the rule's law
    assert 'ocean_depth' in ec.snapshot_tests.__code__.co_names         # the gate reads the run's law, not GDH1
