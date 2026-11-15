"""The crust flow from a measured rheology (planet_tectonics.crust_flow; journal 109)."""
import numpy as np
import pytest

import planet_tectonics.plates as P
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import crust_flow as cf


def _p(rule):
    return {**P.PARAMS, 'crust_flow_rule': rule}


def _kappa(tab, T_km, L_km, grad):
    S = (1-P.PARAMS['rho_crust']/P.PARAMS['rho_mantle'])*grad
    return 10**cf._interp(tab, cf.T_AXIS_KM, cf.L_AXIS_KM, cf.LOGS_AXIS, T_km, L_km, np.log10(S))/grad


def test_geotherm_matches_earths_cratons_and_heats_thick_crust():
    p = P.PARAMS
    qm, q0, tm = cf.moho_state(40e3, 178e3, p)   # a craton: Jaupart & Mareschal 1999, Moho 10-15, surface ~40-45 mW/m^2
    assert 10e-3 <= qm <= 20e-3 and 38e-3 <= q0 <= 55e-3 and 400 < tm < 550
    _, _, t60 = cf.moho_state(60e3, 140e3, p)
    _, _, t40 = cf.moho_state(40e3, 140e3, p)
    assert t60 > t40+200                            # thickening by pure shear heats the base
    _, _, thin = cf.moho_state(40e3, 100e3, p)
    assert thin > t40                               # a thin lithosphere heats it too
    # bulk heat production of a 40 km crust inside Jaupart et al. 2016's 0.65-0.87 uW/m^3
    sh = p['lower_crust_share']
    assert .65e-6 <= p['flow_heat_upper_w_m3']*(1-sh)+p['flow_heat_lower_w_m3']*sh <= .87e-6


@pytest.mark.parametrize('rule', ['channel', 'sheet'])
def test_flux_rises_with_thickness_and_slope_falls_with_lithosphere(rule):
    tab = cf.flux_table(_p(rule))
    assert np.all(np.diff(tab, axis=0) > 0)         # thicker crust: hotter base, more flow
    assert np.all(np.diff(tab, axis=2) > 0)         # steeper surface: more flow
    k = [_kappa(tab, 55, L, .05) for L in (100, 140, 178)]
    assert k[0] > k[1] > k[2]                       # old thick lithosphere keeps the crust cold and stiff
    assert _kappa(tab, 35, 140, .05) < 30           # normal crust over normal lithosphere barely flows


def test_sheet_flows_faster_than_channel():
    a, b = cf.flux_table(_p('channel')), cf.flux_table(_p('sheet'))
    assert np.all(b > a)


def test_edge_flow_conserves_volume():
    g = CubeSphere(24)
    rng = np.random.default_rng(1)
    th = 35e3+20e3*rng.random(g.size)
    dom = rng.random(g.size) < .7
    L = 100.+80.*rng.random(g.size)
    p = _p('sheet')
    K, kmax, capped, counted = cf.edge_kappa(th, dom, L, g.nbr, g.dist, g.axis, cf.flux_table(p), cf.T_AXIS_KM,
                                             cf.L_AXIS_KM, cf.LOGS_AXIS, .14, 30., 2e5)
    assert kmax <= 2e5+30. and counted > 0
    h = .9/(kmax*g.diffusion_bound())
    out = cf.diffuse_edges(th, dom, K, h, g.nbr, g.dist, g.axis, g.area)
    assert np.isclose((out*g.area).sum(), (th*g.area).sum(), rtol=1e-12)
    assert np.array_equal(out[~dom], th[~dom])
    assert out.max() <= th.max() and out[dom].min() >= th[dom].min()   # explicit and stable: no new extremes


@pytest.mark.parametrize('rule', ['channel', 'sheet', 'rey'])
def test_plate_layer_flow_keeps_the_crusts_volume(rule):
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=5, params={'crust_flow_rule': rule, 'threads': 1})
    sim.state()
    for _ in range(3):
        sim.step()
        sim.state()
    w = sim.world
    th = w['thick'].copy()
    th[w['cont']] = th[w['cont']]+np.where(np.arange(w['cont'].sum()) % 7 == 0, 15e3, 0.)   # bumps to spread
    out = cf.flow(sim, th, 1.)
    assert np.isclose((out*g.area).sum(), (th*g.area).sum(), rtol=1e-10)
    assert out[w['cont']].max() < th[w['cont']].max()
    assert sim.erosion_totals['crust_flow_substeps'] >= 1


def test_default_is_rey_and_unknown_rule_raises():
    assert P.PARAMS['crust_flow_rule'] == 'rey'
    sim = P.PlateSim(CubeSphere(16), seed=1, params={'crust_flow_rule': 'nope', 'threads': 1})
    sim.state()
    with pytest.raises(ValueError):
        sim.step()


def test_rey_pins_edges_near_active_convergence_only():
    g = CubeSphere(24)
    rng = np.random.default_rng(2)
    th = 40e3+25e3*rng.random(g.size)
    dom = np.ones(g.size, bool)
    L = np.full(g.size, 140.)
    p = _p('rey')
    fixed, free = cf.flux_table(p, 'channel'), cf.flux_table(p, 'sheet')
    args = (th, dom, L, g.nbr, g.dist, g.axis)
    tail = (cf.T_AXIS_KM, cf.L_AXIS_KM, cf.LOGS_AXIS, .14, 30., 1e30)
    Kf = cf.edge_kappa(*args, fixed, *tail)[0]
    Ks = cf.edge_kappa(*args, free, *tail)[0]
    pin = np.zeros(g.size, bool); pin[:g.size//2] = True
    Kr = cf.edge_kappa(*args, free, *tail, pin, fixed)[0]
    j = np.where(g.nbr >= 0, g.nbr, 0)
    either = pin[:, None] | pin[j]
    ok = g.axis & (g.nbr >= 0)
    assert np.array_equal(Kr[ok & either], Kf[ok & either]) and np.array_equal(Kr[ok & ~either], Ks[ok & ~either])
    assert P.PARAMS['flow_pinned_km'] == 290.


def test_split_step_matches_the_plain_explicit_step():
    g = CubeSphere(48)
    out = {}
    for split in (2e3, 1e30):
        sim = P.PlateSim(g, seed=5, params={'crust_flow_rule': 'sheet', 'threads': 1, 'flow_kappa_split_m2_yr': split})
        sim.state()
        sim.step()
        sim.state()
        w = sim.world
        th = w['thick'].copy()
        th[w['cont']] = th[w['cont']]+np.where(np.arange(w['cont'].sum()) % 5 == 0, 20e3, 0.)
        out[split] = cf.flow(sim, th, 1.)
        assert np.isclose((out[split]*g.area).sum(), (th*g.area).sum(), rtol=1e-10)
    assert np.abs(out[2e3]-out[1e30]).max() < 1e3          # within 1 km of a 20 km bump after 1 Myr
