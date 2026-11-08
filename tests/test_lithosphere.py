"""C, the continental lithosphere's thermal state (fifty-sixth pass, journal 78): the Earth table, the elevation term,
stretching, the plate-layer plumbing, and the C1 review (planet_reviews.coast_budget)."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import lithosphere as lt
from planet_tectonics import plates as P

C1 = DATA/'crust1'


@pytest.mark.skipif(not (C1/'CNtype1-1.txt').exists(), reason='CRUST1.0 not present')
def test_earth_table_and_reference_reproduce_crust1():
    ref, rows = lt.earth_table()
    assert abs(ref-lt.EARTH_LITH_REF_KM) < .01
    assert np.allclose(np.array(rows), np.array(lt.EARTH_LITH_AGE_TABLE))
    assert P.PARAMS['lith_ref_km'] == lt.EARTH_LITH_REF_KM


def test_elevation_term_and_stretching():
    p = P.PARAMS
    t_ref = (p['lith_ref_km']-p['lith_L0_km'])/p['lith_rate_km_myr']
    dz = lt.elevation_term(np.array([t_ref, 0., 3200., np.nan]), p)
    assert abs(dz[0]) < 1e-9 and dz[3] == 0.
    assert abs(dz[1]-30/3300*(130.67-93.6)*1e3) < 1e-6                 # young lithosphere stands ~337 m higher
    assert abs(dz[2]+30/3300*(93.6+128.-130.67)*1e3) < 1e-6            # an early Archean root ~827 m lower
    a = np.array([1000., 1000., 50.])
    s = lt.stretched_age(a, np.array([1., .5, .5]), p)
    assert abs(s[0]-1000.) < 1e-9
    assert s[1] == 0. and s[2] == 0.                                   # L 133.6 -> 66.8 km: below L0, age 0
    assert abs(lt.stretched_age(np.array([3000.]), np.array([.8]), p)[0]-((93.6+120.)*.8-93.6)/.04) < 1e-9


def test_start_ages_follow_distance_and_craton():
    d = np.linspace(0., 3000e3, 1000)
    young = lt.start_ages(d, np.zeros(1000))
    old = lt.start_ages(d, np.ones(1000))
    mid = lt.start_ages(d, np.full(1000, .5))
    assert young.min() == 33. and old.max() == 3200.
    assert np.all(old >= young)
    assert mid[d >= 1000e3].min() > mid[d < 50e3].max()                 # young rims, old cores
    tab = np.array(lt.EARTH_LITH_AGE_TABLE)
    assert np.all(np.diff(tab, axis=1) >= 0)                           # quantiles rise within a row


def test_plate_layer_carries_and_resets_thermal_ages():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'lithosphere': 'thermal_age'})
    w = sim.world
    assert not np.isnan(w['therm'][w['cont']]).any()
    assert np.isnan(w['therm'][~w['cont']]).all()
    for _ in range(25):
        sim.step()
    st = sim.state()
    c = st['continental']
    la = st['lithosphere_age_myr']
    assert not np.isnan(la[c]).any()                                   # every continental node carries an age
    assert np.isnan(la[~c]).all()
    assert (la[c] < 25.).any()                                         # resets happened
    assert (la[c] > 1000.).any()                                       # old interiors kept theirs (aged by 25 Myr)


def test_elevation_term_is_local_and_off_with_none():
    g = CubeSphere(48)                                                 # 'thermal_age' is PARAMS since pass 69
    sim = P.PlateSim(g, seed=3, params={'erosion': False, 'threads': 1, 'lithosphere': 'none'})   # exact compare on one thread (pass 57:
    sim.elevation(); sp = sim.sea_level                    # all cores' flexure differs at ~1e-15 m between calls; the
                                                           # loads read the previous sea level)
    z0 = sim.elevation(); s0 = sim.sea_level
    sim.world['therm'] = np.where(sim.world['cont'], sim.t-3000., np.nan)
    sim.sea_level = sp
    assert np.array_equal(sim.elevation(), z0)                         # 'none': the age does nothing
    sim.p['lithosphere'] = 'thermal_age'
    sim.sea_level = sp
    z1 = sim.elevation(); s1 = sim.sea_level
    dl = lt.elevation_term(np.full(g.size, 3000.), sim.p)
    dry = sim.world["cont"] & (z0+sp-s0 > 1200.)                      # still above the sea with the term
    assert np.allclose((z1+s1-(z0+s0))[dry], dl[dry])                # dry land moves by the term exactly
    assert np.allclose((z1+s1-(z0+s0))[~sim.world['cont']], 0.)        # ocean floor untouched (on the datum)


def test_coast_budget_reads_a_fall_below_the_airy_line():
    from planet_reviews import coast_budget as cb
    g = CubeSphere(217)
    ang = np.arccos(np.clip(g.xyz @ np.array([0., 0., 1.]), -1, 1))
    land = ang < np.radians(25.)
    R = g.radius/1e3
    d = (np.radians(25.)-ang)*R                                        # km inland (approximately)
    p = P.PARAMS
    Tc = np.where(land, 36e3+5.*d, 7e3); S = np.zeros(g.size)
    on = cb.airy_m(Tc, S, p)
    z = np.where(land, on, -4000.)
    b = cb.budget(g, z, Tc, S, land, p)
    assert abs(b['c1_m']) < 1e-6 and b['rise_m'] > 0                    # on its Airy line: no fall
    z2 = np.where(land, on-np.clip(d, 0, 900)*.8, -4000.)              # falls 0.8 m per km inland below it
    b2 = cb.budget(g, z2, Tc, S, land, p)
    assert -800. < b2['c1_m'] < -500.


def test_coast_budget_keep_and_term():
    """keep restricts the nodes read (distances still from the whole land mask); term gives the band means' change."""
    from planet_reviews import coast_budget as cb
    g = CubeSphere(217)
    ang = np.arccos(np.clip(g.xyz @ np.array([0., 0., 1.]), -1, 1))
    land = ang < np.radians(25.)
    d = (np.radians(25.)-ang)*g.radius/1e3
    p = P.PARAMS
    Tc = np.where(land, 38e3, 7e3); S = np.zeros(g.size)
    z = np.where(land, cb.airy_m(Tc, S, p)+np.where(d < 50, 300., 0.), -4000.)    # the coast band stands 300 m high
    term = np.where(land, -.2*np.clip(d, 0, 1200), 0.)
    b = cb.budget(g, z, Tc, S, land, p, term=term)
    assert abs(b['c1_m']+300.) < 1e-6 and -250. < b['term_fall_m'] < -150.
    keep = d >= 25.                                                    # drop the outer half of the coast band
    bk = cb.budget(g, z, Tc, S, land, p, keep=keep)
    assert bk['bins'][0]['area_share'] < b['bins'][0]['area_share'] and abs(bk['c1_m']+300.) < 1e-6


def test_elevation_gate_reads_c_in_the_columns_own_level(tmp_path):
    """The thin-crust gate (elevation_check) checks the column against the model's own column, C's term included
    (pass 69): a C run's state reads as the same state without C, and the old reading (C left out) does not."""
    from planet_reviews import elevation_check as ec
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'erosion': False, 'threads': 1, 'lithosphere': 'thermal_age'})
    sim.state()
    for _ in range(3):
        sim.step(); sim.state()
    np.savez(tmp_path/'c.npz', **sim.state())
    pc = dict(sim.p)
    sim.p['lithosphere'] = 'none'
    np.savez(tmp_path/'n.npz', **sim.state())
    pn = dict(sim.p)
    c, n = ec.snapshot_tests(tmp_path/'c.npz', g, pc), ec.snapshot_tests(tmp_path/'n.npz', g, pn)
    old = ec.snapshot_tests(tmp_path/'c.npz', g, pn)                   # C's state read without C's term
    assert c['thin_crust_area'] > 0
    assert abs(c['thin_crust_offset_m']-n['thin_crust_offset_m']) < 30.
    assert old['thin_crust_offset_m']-n['thin_crust_offset_m'] > 100.
    assert c['young_ocean_offset_m'] == old['young_ocean_offset_m']   # the ocean side does not read C
