"""The reviews' label-free continent (twenty-second pass, agreed 25 Oct): total thickness >= 20 km."""
import json

import numpy as np
import pytest

from planet_reviews import crust_review as cr
from planet_grid.paths import DATA

C1 = DATA/'crust1'


@pytest.mark.skipif(not (C1/'crsthk.xyz').exists(), reason='CRUST1.0 not present')
def test_label_free_rule_reproduces_crust1_types_on_earth():
    e = json.loads((C1/'crust1-continental-summary.json').read_text(encoding='utf-8'))
    m = cr.earth_label_free()
    assert abs(m['area']-e['continental_area_fraction']) < .005
    assert abs(m['thin_lt30']-e['continental_share_thinner_than_30km']) < .005
    assert abs(m['thick_gt45']-e['continental_share_thicker_than_45km']) < .005
    assert abs(m['median_km']-e['continental_total_thickness_km_percentiles']['50']) <= .5


def test_column_thickness_counts_sediment_on_both_sides():
    z = {'continental': np.array([1, 1, 0, 0], np.uint8),
         'crust_thickness_m': np.array([35e3, 18e3, 7e3, 15e3]), 'sediment_m': np.array([2e3, 3e3, 1e3, 6e3])}
    with_layer = cr.column_thickness_m(z, {'continental_sediment': True})
    assert np.allclose(with_layer, [37e3, 21e3, 8e3, 21e3])
    # before the sediment layer, deposits on continents were booked as crust: the sediment field there is not a layer
    assert np.allclose(cr.column_thickness_m(z, {'continental_sediment': False}), [35e3, 18e3, 8e3, 21e3])
    # a thinned margin (flag continental, 18 km) is ocean to the rule; a thick sediment pile on ocean crust is continent
    assert list(cr.label_free_continent(with_layer)) == [True, True, False, True]


def test_coast_profile_measures_distance_inland_and_the_dome():
    """coast_profile (pass 54): on a spherical-cap continent of angular radius r the deepest node lies ~r inland
    (great circle, km); a dome (elevation rising inland) reads a positive rise and no far lowland, a flat 100 m
    continent reads zero rise and all far land below 200 m."""
    from planet_grid.sphere_graph import CubeSphere
    from planet_reviews import coast_profile as cp
    g = CubeSphere(217)                                                   # the review's grid (0-50 km bin filled)
    ang = np.arccos(np.clip(g.xyz @ np.array([0., 0., 1.]), -1, 1))
    r = np.radians(30.)
    land = ang < r
    d = cp.distance_to_ocean_km(g, land)
    R = g.radius/1e3
    assert d[~land].max() == 0.
    assert abs(d[land].max() - r*R) < 2.*np.sqrt(g.area.mean())/1e3      # within ~2 grid spacings
    dome = np.where(land, 2.*(r-ang)*R, -4000.)                           # 2 m per km inland
    p = cp.profile(g, dome, land)
    assert p['rise_m'] > 1000. and p['lowland_far'] == 0.
    flat = np.where(land, 100., -4000.)
    q = cp.profile(g, flat, land)
    assert q['rise_m'] == 0. and q['lowland_far'] == 1.
    assert abs(sum(b['area_share'] for b in q['bins']) - 1.) < 1e-12
