"""flow_power (ninety-ninth pass, journal 132): D's slab-flow pattern at Earth's measured power per degree."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_tectonics import plates as P
from planet_tectonics.mantle import measured_power_air_m2

pytestmark = pytest.mark.skipif(not (DATA/'mantle'/'visc.sh08').exists(), reason='data/mantle/visc.sh08 missing')

FLOW_P = {'rho_mantle': 3300., 'flow_lmax': 31, 'flow_z_cut_km': 300., 'slab_dip_deg': 58., 'lm_sink_mm_yr': 12.,
          'viscosity_file': None, 'threads': 1, 'slab_alpha': 'thermo'}


def _fill(fl, g, seed=2):
    """A few trench segments' parcels, advanced until most lie below the cut."""
    rng = np.random.default_rng(seed)
    for k in range(6):
        I = rng.choice(g.size, 40, replace=False)
        d = np.cross(g.xyz[I], rng.normal(size=3)); d /= np.linalg.norm(d, axis=1)[:, None]
        fl.record.deposit(g.xyz[I], d, np.full(40, 1e16), np.full(40, .05), -100.+k)
        fl.record.advance(5., g.radius)
    fl.record.advance(30., g.radius)
    fl.record.merge(g, fl.z_cut)


def _pair():
    g = CubeSphere(24)
    a = mf.MantleFlow(g, {**FLOW_P, 'flow_power': 'flow'})
    b = mf.MantleFlow(g, {**FLOW_P, 'flow_power': 'measured'})
    _fill(a, g); _fill(b, g)
    a.update(); b.update()
    return g, a, b


def test_measured_power_takes_earths_power_per_degree_and_keeps_the_pattern():
    g, a, b = _pair()
    lm = a.lmax
    assert np.allclose(a.raw_power, b.raw_power, rtol=1e-12)          # the same record, the same flow
    pw = mf.spectrum(b.H, lm)*4*np.pi
    target = measured_power_air_m2()
    assert (a.raw_power[1:31] > 0).all()
    assert np.allclose(pw[1:31], target[1:31], rtol=1e-9)              # Earth's power per degree, l 1..30
    assert pw[0] == 0. and pw[31] == 0.                                # no degree 0; none beyond the measurement
    L = mf.alm_l(lm)
    for l in (1, 2, 5, 17, 30):                                        # the shape within each degree is the flow's
        ha, hb = a.H[L == l], b.H[L == l]
        assert np.allclose(hb/np.abs(hb).max(), ha/np.abs(ha).max(), atol=1e-12)
    rms = np.sqrt(np.mean(b.h_air**2))                                 # 0.41 km air-loaded (mantle.py's header)
    assert 300. < rms < 520.
    assert np.allclose(a.V, b.V)                                       # the surface velocity is the flow's own


def test_lows_stay_over_the_slabs():
    g, a, b = _pair()
    rec = b.record
    deep = rec.z >= b.z_cut
    near = np.zeros(g.size)
    from planet_grid import fields
    np.add.at(near, fields.nearest_nodes(g, rec.pos[deep]), rec.mass[deep])
    # the field under the record's mass is below its mean in both, and the two fields are positively correlated
    w = near/near.sum()
    assert (w*a.h_air).sum() < a.h_air.mean() and (w*b.h_air).sum() < b.h_air.mean()
    assert np.corrcoef(a.h_air, b.h_air)[0, 1] > .3


def test_an_empty_record_gives_no_field_and_bad_values_raise():
    g = CubeSphere(24)
    b = mf.MantleFlow(g, {**FLOW_P, 'flow_power': 'measured'})
    b.update()
    assert not np.any(b.H) and not np.any(b.h_air)
    with pytest.raises(ValueError):
        mf.MantleFlow(g, {**FLOW_P, 'flow_power': 'scaled'})


def test_default_is_the_placed_field_and_a_measured_run_saves_the_raw_power():
    # adopted 21 Nov (the hundredth pass, journal 134: the Placed Field Sweep)
    assert P.PARAMS['flow_power'] == 'measured' and P.PARAMS['dyn_rule'] == 'slabs'
    g = CubeSphere(24)
    s = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'erosion': False, 'threads': 1, 'dyn_rule': 'slabs',
                                      'flow_power': 'measured'})
    for _ in range(3):
        s.step()
    st = s.state()
    assert st['flow_power_raw_m2'].shape == (s.mantle.flow.lmax+1,)
    d = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'erosion': False, 'threads': 1, 'dyn_rule': 'slabs',
                                      'flow_power': 'flow'})
    d.step()
    assert 'flow_power_raw_m2' not in d.state()


def test_the_start_is_the_controls_and_hands_over_to_the_placed_field():
    """The record starts empty: the start's mantle is dyn_rule 'measured''s field (the same stream), so the start is the
    control's; it hands over to the placed flow at constant variance within one knot interval (dt_knot_myr)."""
    g = CubeSphere(24)
    kw = {'water_volume_km3': None, 'threads': 1}
    a = P.PlateSim(g, seed=3, params={**kw, 'dyn_rule': 'measured'})
    b = P.PlateSim(g, seed=3, params={**kw, 'dyn_rule': 'slabs', 'flow_power': 'measured'})
    sa, sb = a.state(), b.state()
    for k in ('elevation_m', 'sediment_m', 'crust_thickness_m', 'dynamic_topography_m', 'dynamic_topography_air_m'):
        assert np.array_equal(sa[k], sb[k]), k
    m = b.mantle
    t1 = m.t0+P.PARAMS['dt_knot_myr']
    assert np.array_equal(m.air(t1), m.flow.h_air) and np.array_equal(m.air(m.t0), m._measured_at(m.t0))
    with pytest.raises(ValueError):
        from planet_tectonics.mantle import Mantle
        Mantle(g, np.random.default_rng(1), -100., params={'dyn_rule': 'slabs', 'flow_power': 'measured'})
