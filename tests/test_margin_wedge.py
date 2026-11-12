"""S, the margin wedge (sixty-first pass, journal 83): what overflows the graded shelf is deposited in front of its
source up to Earth's measured margin profile before any rest fills the sea-floor lows. Build checks: the wedge fills
from the source outward and every m^3 is placed; a supply larger than the wedge holds leaves the rest to the lows as
before; off, nothing changes; the plate layer's ceiling applies only where the profile was measured."""
import numpy as np

from planet_grid import sphere_graph as sg
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import plates as P


def _strip():
    g = sg.flat(1, 12, 10e3)
    z = np.full(12, -5000.); z[0] = 100.
    ocean = z < 0
    sea_in = np.zeros(12)
    d = gl.distance_to(~ocean, g.nbr, g.dist, 15e3)
    top = np.full(12, -1e9); top[2:] = np.linspace(-2000., -4500., 10)
    return g, z, ocean, sea_in, d, top


def _run(g, z, ocean, sea_in, d, top):
    return gl.marine_sediment_kept(z, ocean, sea_in, g.nbr, g.dist, g.area, 10., 1.2e-3, d, 15e3, top)


def test_wedge_fills_from_the_source_outward_and_places_every_m3():
    g, z, ocean, sea_in, d, top = _strip()
    delta = (2.-z[1])*g.area[1]                            # the delta node at the river mouth fills to +2 m first
    cap = ((top[2:]-z[2:])*g.area[2:]).sum()
    sea_in[1] = delta+.5*cap
    dep, left = _run(g, z, ocean, sea_in, d, top)
    assert np.isclose((dep*g.area).sum(), sea_in[1]) and left.sum() == 0.
    surf = z+dep
    full = np.isclose(surf[2:], top[2:])
    k = np.argmin(full)                                    # the fill front
    assert full[:k].all() and k >= 2                       # contiguous from the source outward
    assert np.allclose(dep[2+k+1:], 0.)                    # nothing beyond the front
    assert (surf[2:] <= top[2:]+1e-9).all()                # never above the profile


def test_supply_beyond_the_wedge_goes_to_the_lows_as_before():
    g, z, ocean, sea_in, d, top = _strip()
    cap = (2.-z[1])*g.area[1]+((top[2:]-z[2:])*g.area[2:]).sum()
    sea_in[1] = cap+3e11
    dep, left = _run(g, z, ocean, sea_in, d, top)
    assert np.allclose(z[2:]+dep[2:], top[2:])             # the wedge is full to the profile
    assert np.isclose(left.sum(), 3e11) and np.isclose((dep*g.area).sum()+left.sum(), sea_in[1])
    z1 = gl.fill_and_spill(z+dep, left, g.nbr, g.area)
    assert np.isclose(((z1-z-dep)*g.area).sum(), 3e11)     # the rest is placed by fill_and_spill


def test_turbidite_path_builds_the_wedge_from_a_pit():
    """fill_and_spill with the ceiling: a remainder left in a shallow pit fills the pit to its sill, spills, and
    builds the wedge down the slope before ponding in the deep low; every m^3 is placed, off is the old rule."""
    g = sg.flat(1, 9, 10e3)
    z = np.array([-400., -500., -300., -2500., -3000., -3500., -4000., -4500., -4400.])
    ceil = np.array([-1e9, -1e9, -1e9, -2000., -2500., -3000., -3500., -4000., -1e9])
    vol = np.zeros(9); vol[1] = 1e11                       # 1 km over a 1e8 m^2 cell
    z0 = gl.fill_and_spill(z, vol, g.nbr, g.area)
    assert np.array_equal(z0, gl.fill_and_spill(z, vol, g.nbr, g.area, np.full(9, -1e9)))
    z1, wv = gl.fill_and_spill(z, vol, g.nbr, g.area, ceil, return_wedge=True)
    assert np.isclose(((z1-z)*g.area).sum(), 1e11) and np.isclose(((z0-z)*g.area).sum(), 1e11)
    assert np.allclose(z1[:3], -300.)                      # the pit fills to its sill (3e10 m^3) and spills
    assert np.isclose(z1[3], -2000.) and np.isclose(z1[4], -2800.)   # the rest builds the slope from the top down
    assert np.isclose(wv.sum(), 7e10) and np.allclose(z1[5:], z[5:])
    assert z0[3] == z[3]                                   # without it the slope is bypassed


def test_off_is_the_old_rule_exactly():
    g, z, ocean, sea_in, d, top = _strip()
    sea_in[1] = 1e13
    a = _run(g, z, ocean, sea_in, d, None)
    b = _run(g, z, ocean, sea_in, d, np.full(12, -1e9))
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
    assert np.allclose(a[0][2:], 0.)                       # without the wedge nothing stays beyond the shelf


def test_plate_layer_ceiling_only_where_measured_and_ledger_closes():
    g = CubeSphere(96)                                     # at n = 48 (192 km) the young floor stands above the profile
    # crust_density_rule and lith_frame pinned to the trajectory this test was set on (its seed grows no wedge with
    # 'pressure', nor booked spill with 'lamb')
    sim = P.PlateSim(g, seed=3, params={'margin_wedge': True, 'threads': 1, 'crust_density_rule': 'uniform',
                                        'lith_frame': 'artemieva'})
    for _ in range(12):
        sim.step()
    p, w = sim.p, sim.world
    top = sim._wedge_top(w['thick'])
    on = top > -1e8
    assert on.any()
    from planet_tectonics.isostasy import sediment_column
    col = w['thick']+sediment_column(w['sed'], w['cont'], p)[0]
    assert not (on & (col >= 20e3)).any()                  # never on the label-free continent (thinned margins may)
    assert top[on].max() <= -2083.+1e-6 and top[on].min() >= -4761.-1e-6
    near = gl.distance_to(sim.btype != 0, g.nbr, g.dist, 300e3) <= 300e3
    assert not (on & near).any()                           # not within 300 km of a plate boundary
    z = sim.elevation(w['thick'], w['sed'])
    lp = {**gl.DEFAULTS, 'open_edges': False}
    K = np.full(g.size, lp['K'])
    r = gl.surface_step(g, z, z < -1500., None, np.ones(g.size, bool), K, lp, 2.5e5, sim.h_min, sim.bound,
                        keep_deep=True, wedge_top=top)
    r0 = gl.surface_step(g, z, z < -1500., None, np.ones(g.size, bool), K, lp, 2.5e5, sim.h_min, sim.bound,
                         keep_deep=True)
    gap = ((r['z']-z)*g.area).sum()+r['deep_m3'].sum()
    gap0 = ((r0['z']-z)*g.area).sum()+r0['deep_m3'].sum()   # the kernel's own residual (~4e-6 of the eroded volume
    assert r['wedge_m3'].sum() > 0.                          # here, journal 83), not the wedge's
    assert abs(gap-gap0) < 1e-9*(r['eroded']*g.area).sum()   # every m^3 the wedge takes it takes from the lows
    # the lows lose at least what the wedge holds (more where a wedge-filled path turns back to a deep node inside
    # the 150 km shelf zone, which at coarse spacing can lie far below the shelf's top)
    assert r0['deep_m3'].sum()-r['deep_m3'].sum() >= r['wedge_m3'].sum()*(1-1e-9)
    sim._erode(w['thick'].copy())
    t = sim.erosion_totals
    assert t['margin_wedge_m3'] > 0. and t['margin_wedge_spill_m3'] >= 0.   # booked in the plate layer's totals


def test_default_is_off():
    assert P.PARAMS['margin_wedge'] is False
    km, depth = np.array(P.PARAMS['margin_wedge_profile']).T
    assert (np.diff(km) > 0).all() and (np.diff(depth) < 0).all()   # the rise deepens seaward
