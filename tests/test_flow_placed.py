"""flow_placed (pass 101, journal 136): the slab record places only the share of each degree's measured power that the
same flow on Earth's own subduction history accounts for on Earth (mantle.EARTH_SLAB_CORR); the rest is plate-blind."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_tectonics import plates as P
from planet_tectonics.mantle import EARTH_SLAB_CORR, Mantle, measured_power_air_m2, placed_share

pytestmark = pytest.mark.skipif(not (DATA/'mantle'/'visc.sh08').exists(), reason='data/mantle/visc.sh08 missing')

FLOW_P = {'rho_mantle': 3300., 'flow_lmax': 31, 'flow_z_cut_km': 300., 'slab_dip_deg': 58., 'lm_sink_mm_yr': 12.,
          'viscosity_file': None, 'threads': 1, 'slab_alpha': 'thermo', 'flow_power': 'measured'}


def _fill(fl, g, seed=2):
    rng = np.random.default_rng(seed)
    for k in range(6):
        I = rng.choice(g.size, 40, replace=False)
        d = np.cross(g.xyz[I], rng.normal(size=3)); d /= np.linalg.norm(d, axis=1)[:, None]
        fl.record.deposit(g.xyz[I], d, np.full(40, 1e16), np.full(40, .05), -100.+k)
        fl.record.advance(5., g.radius)
    fl.record.advance(30., g.radius)
    fl.record.merge(g, fl.z_cut)


def test_the_share_is_earths_measured_correlation_squared():
    sh = placed_share(31)
    assert len(EARTH_SLAB_CORR) == 30 and sh.shape == (32,)
    assert sh[0] == 0. and sh[31] == 0.
    r = np.asarray(EARTH_SLAB_CORR)
    assert np.allclose(sh[1:31], np.clip(r, 0, None)**2)
    assert sh[1] > .98 and .7 < sh[2] < .8 and sh[3] < .15                 # l 1-2 placed; little of l 3
    pw = measured_power_air_m2()
    assert .15 < (pw*sh[:31]).sum()/pw.sum() < .3                          # the placed share of Earth's whole power
    assert (pw[4:]*sh[4:31]).sum()/pw[4:].sum() < .05                      # degrees 4-30: hardly any


def test_the_flow_carries_the_share_of_the_measured_power_with_the_same_pattern():
    g = CubeSphere(24)
    a = mf.MantleFlow(g, FLOW_P)
    b = mf.MantleFlow(g, {**FLOW_P, 'flow_placed': 'earth_corr'})
    _fill(a, g); _fill(b, g)
    a.update(); b.update()
    lm = a.lmax
    pa, pb = mf.spectrum(a.H, lm)*4*np.pi, mf.spectrum(b.H, lm)*4*np.pi
    sh = placed_share(lm)
    assert np.allclose(pb[1:31], pa[1:31]*sh[1:31], rtol=1e-9, atol=1e-6)
    L = mf.alm_l(lm)
    for l in (1, 2, 3):                                                    # the shape within each degree is the flow's
        ha, hb = a.H[L == l], b.H[L == l]
        assert np.allclose(hb/np.abs(hb).max(), ha/np.abs(ha).max(), atol=1e-12)
    assert not np.any(b.H[L == 4])                                         # r_4 < 0: nothing placed
    assert np.allclose(a.V, b.V)


def test_the_start_is_the_controls_and_the_rest_is_plate_blind():
    g = CubeSphere(24)
    kw = {'water_volume_km3': None, 'threads': 1, 'dyn_rule': 'slabs', 'flow_power': 'measured'}
    a = P.PlateSim(g, seed=3, params={**kw, 'flow_placed': 'all'})
    b = P.PlateSim(g, seed=3, params={**kw, 'flow_placed': 'earth_corr'})
    sa, sb = a.state(), b.state()
    for k in ('elevation_m', 'sediment_m', 'crust_thickness_m', 'dynamic_topography_m', 'dynamic_topography_air_m'):
        assert np.array_equal(sa[k], sb[k]), k
    m = b.mantle
    assert a.mantle.unplaced is None and len(m.unplaced) == len(m.measured)
    for x, y in zip(a.mantle.measured, m.measured):                        # the start's field: the same stream
        assert np.array_equal(x, y)
    t1 = m.t0+P.PARAMS['dt_knot_myr']
    assert np.allclose(m.air(t1), m.flow.h_air+m._knots_at(m.unplaced, t1))
    pw = measured_power_air_m2()
    expect = np.sqrt((pw*(1.-placed_share(pw.size-1))).sum()/(4*np.pi))
    u = np.concatenate(m.unplaced)
    assert .7*expect < np.sqrt(np.mean(u**2)) < 1.3*expect                 # the power the slabs leave
    assert abs(np.corrcoef(m.unplaced[0], m.measured[0])[0, 1]) < .3       # independent of the start's field


def test_default_and_bad_values():
    assert P.PARAMS['flow_placed'] == 'all'
    g = CubeSphere(24)
    with pytest.raises(ValueError):
        Mantle(g, np.random.default_rng(1), -100., params={'dyn_rule': 'slabs', 'flow_power': 'flow',
                                                           'flow_placed': 'earth_corr'}, field_rng=np.random.default_rng(2))
    with pytest.raises(ValueError):
        Mantle(g, np.random.default_rng(1), -100., params={'dyn_rule': 'slabs', 'flow_power': 'measured',
                                                           'flow_placed': 'some'}, field_rng=np.random.default_rng(2))
