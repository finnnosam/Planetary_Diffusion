"""rift_traction (pass 103, journal 139): the split rule's forces read the basal traction of the planet's own mantle flow,
the flow scaled with its topography (MantleFlow.placed_velocity)."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_tectonics import plates as P
from planet_tectonics.motion import YR_S

from test_flow_power import FLOW_P, _fill
from test_rift_line import _big_continental_plate, _sim


def test_default_on():
    assert P.PARAMS['rift_traction'] is True                 # default since pass 104 (journal 140)


@pytest.mark.skipif(not (DATA/'mantle'/'visc.sh08').exists(), reason='data/mantle/visc.sh08 missing')
def test_placed_velocity_is_the_scaled_flow():
    """'flow': the flow's own velocity; 'measured': every degree's velocity takes its topography's factor, so the
    velocity and the topography of one degree scale together (the record with that degree's loads scaled)."""
    g = CubeSphere(24)
    a = mf.MantleFlow(g, {**FLOW_P, 'flow_power': 'flow'})
    b = mf.MantleFlow(g, {**FLOW_P, 'flow_power': 'measured'})
    _fill(a, g); _fill(b, g)
    a.update(); b.update()
    assert np.array_equal(a.placed_velocity(), a.v_xyz)
    L = mf.alm_l(b.lmax)
    ok = np.abs(a.H) > 0
    ratio = np.abs(b.H[ok])/np.abs(a.H[ok])                         # the topography's factor per coefficient
    assert np.allclose(ratio, b.scale[L[ok]], rtol=1e-9)
    vs = b.placed_velocity()
    for l in (1, 2, 7):                                             # one degree alone: velocity x its factor
        one = np.where(L == l, b.V, 0)
        v1 = mf.gradient(one*b.scale[l], b.lmax, g.xyz)
        v0 = mf.gradient(one, b.lmax, g.xyz)
        assert np.allclose(v1, b.scale[l]*v0, rtol=1e-9, atol=1e-15)
    assert np.allclose(vs, mf.gradient(b.V*b.scale[L], b.lmax, g.xyz))
    assert np.allclose((vs*g.xyz).sum(1), 0., atol=1e-12)           # tangential


class _Flow:
    """A stand-in flow: radially divergent from the point c at speed u (m/yr), zero at c and its antipode."""
    def __init__(self, g, c, u):
        t = c[None, :]-(g.xyz@c)[:, None]*g.xyz                     # toward c, tangential
        self.v = -u*t/np.maximum(np.linalg.norm(t, axis=1), 1e-12)[:, None]
        self.v *= np.sin(np.arccos(np.clip(g.xyz@c, -1, 1)))[:, None]

    def placed_velocity(self):
        return self.v


def _with_flow(traction, u=.05):
    sim, pl, own = _big_continental_plate('force')
    g = sim.g
    c = (g.area[own, None]*g.xyz[own]).sum(0); c /= np.linalg.norm(c)
    sim.mantle.flow = _Flow(g, c, u)
    sim.p['rift_traction'] = traction
    sim._cut_forces_t = None
    return sim, pl, own


def test_traction_adds_the_basal_drag_on_the_flow():
    s0, pl, own = _with_flow(False)
    s1, _, _ = _with_flow(True)
    t0, w0 = s0._cut_forces()[:2]
    t1, w1 = s1._cut_forces()[:2]
    g, p = s1.g, s1.p
    want = np.cross(g.xyz, s1.mantle.flow.placed_velocity()/YR_S)*(p['basal_drag_pa_s_m']*g.radius*w1)[:, None]
    assert np.allclose(t1-t0, want, rtol=1e-12, atol=0)
    assert np.array_equal(w0, w1)


def test_divergent_flow_opens_the_overriding_supercontinent():
    """The test supercontinent overrides all round: without traction no cut is in tension (the cheapest is taken,
    test_rift_line); over a flow diverging from its centre the force cuts are in tension and break."""
    kinds = {}
    for traction in (False, True):
        sim, pl, own = _with_flow(traction)
        sim.rng = np.random.default_rng(11)
        sim.chord_stats = {'chord': 0, 'bisector': 0}
        for _ in range(8):
            lab = sim._rift_line(pl, own)
            if sim._line_kind == 'force':
                assert sim._cut_tension(own, own & (lab == 1)) > 0
        kinds[traction] = dict(sim.chord_stats)
    assert kinds[False].get('force', 0) == 0, kinds
    assert kinds[True].get('force', 0) >= 6, kinds


def test_without_a_flow_traction_is_moot():
    sim = _sim({'rift_line': 'force', 'lithosphere': 'thermal_age', 'dyn_rule': 'measured'})
    assert sim.mantle.flow is None                                  # no slab flow: nothing to add
    t0 = sim._cut_forces()[0].copy()
    sim.p['rift_traction'] = True; sim._cut_forces_t = None
    assert np.array_equal(sim._cut_forces()[0], t0)
