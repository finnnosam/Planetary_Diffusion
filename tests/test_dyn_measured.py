"""dyn_rule 'measured' (fifty-eighth pass, journal 80): the dynamic topography from Earth's measured spectrum
(planet_tectonics.mantle.VD2020_POWER_KM2, Valentine & Davies 2020) and fields.spectrum_field."""
import numpy as np
import pytest

from planet_grid import fields
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle as M
from planet_tectonics import plates as P

VD = DATA/'mantle'/'vd2020_sph_coeff_high_accuracy_spot.dat'


def test_m1_the_table_is_the_published_coefficients():
    """VD2020_POWER_KM2 = sum over m of (mean^2 + standard error^2) of the published coefficients, l = 1 .. 30."""
    if not VD.exists():
        pytest.skip('data/mantle/vd2020_sph_coeff_high_accuracy_spot.dat missing')
    d = np.loadtxt(VD)
    l = d[:, 0].astype(int)
    p = np.array([(d[l == k, 2]**2+d[l == k, 3]**2).sum() for k in range(1, 31)])
    assert l.max() == 30 and len(l) == 31**2-1
    assert np.allclose(p, M.VD2020_POWER_KM2, atol=6e-5)
    # the published statements: peak at l 2 (likely 0.46-0.76 km^2), ~0.02 at l 30; ~85% in the first 15 degrees
    assert np.argmax(p) == 1 and .46 <= p[1] <= .76 and abs(p[29]-.02) < .005
    assert abs(p[:15].sum()/p.sum()-.86) < .01 and abs(p[:3].sum()/p.sum()-.287) < .005


def test_m2_spectrum_field_has_the_requested_power():
    g = CubeSphere(48)
    w = g.area/g.area.sum()
    pw = np.zeros(31); pw[5] = 4*np.pi            # one degree, mean square 1
    rng = np.random.default_rng(0)
    ms = [float((w*fields.spectrum_field(g.xyz, pw, rng)**2).sum()) for _ in range(200)]
    assert abs(np.mean(ms)-1) < .05               # chi-square with 11 dof: the mean of 200 is within ~2%
    pw = M.measured_power_air_m2()
    ms = [float((w*fields.spectrum_field(g.xyz, pw, rng)**2).sum()) for _ in range(200)]
    want = pw.sum()/(4*np.pi)                      # (405 m)^2
    assert abs(np.sqrt(want)-405.4) < 1 and abs(np.mean(ms)/want-1) < .06
    f = fields.spectrum_field(g.xyz, pw, rng)
    assert abs((w*f).sum()) < 1e-3*np.sqrt(want)   # no degree 0 (to the cube cells' quadrature)


def test_m3_constant_variance_blend_and_water_loading():
    g = CubeSphere(24)
    rng, frng = np.random.default_rng([3, 7]), np.random.default_rng([3, 8])
    m = M.Mantle(g, rng, -100., params={'dyn_rule': 'measured', 'rho_mantle': 3300., 'rho_water': 1030.}, field_rng=frng)
    k0, k1 = m.knots_t[1], m.knots_t[2]
    assert np.allclose(m.air(k0), m.measured[1]) and np.allclose(m.air(k1), m.measured[2])
    assert np.allclose(m.air(.5*(k0+k1)), (m.measured[1]+m.measured[2])/np.sqrt(2))
    wet = np.zeros(g.size, bool); wet[::2] = True
    t = k0+7.
    h = m.topography(t, wet)-m.swells(t)
    assert np.allclose(h[wet], m.air(t)[wet]*3300/2270) and np.allclose(h[~wet], m.air(t)[~wet])
    # the measured fields come from their own stream: the hotspots are those of the random rule
    r = M.Mantle(g, np.random.default_rng([3, 7]), -100.)
    assert [h['pos'].tolist() for h in r.hotspots] == [h['pos'].tolist() for h in m.hotspots]
    with pytest.raises(ValueError):
        M.Mantle(g, np.random.default_rng(1), -100., params={'dyn_rule': 'measured'})


def test_m4_plate_sim_runs_with_the_measured_field_now_the_default():
    assert P.PARAMS['dyn_rule'] == 'measured'           # adopted 6 Nov (journal 88)
    g = CubeSphere(24)
    a = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'erosion': False, 'threads': 1, 'dyn_rule': 'random'})
    b = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'threads': 1, 'dyn_rule': 'measured'})
    assert a.mantle.measured is None and b.mantle.measured is not None and b.mantle.flow is None
    assert [h['pos'].tolist() for h in a.mantle.hotspots] == [h['pos'].tolist() for h in b.mantle.hotspots]
    for _ in range(3):
        a.step(); b.step()
    sa, sb = a.state(), b.state()
    assert 'dynamic_topography_air_m' in sb and 'dynamic_topography_air_m' not in sa
    # with erosion on, the last step's surface sets the water factor (z_last < 0: x 3300 / 2270, swells unloaded)
    sw = b.mantle.swells(b.t)
    r = (sb['dynamic_topography_m']-sw)/(sb['dynamic_topography_air_m']-sw)
    ok = np.abs(sb['dynamic_topography_air_m']-sw) > 1.
    assert np.allclose(np.unique(np.round(r[ok], 3)), [1., round(3300/2270, 3)], atol=2e-3)
