"""Pass 76 (journal 100): crust_density_rule 'pressure', the lower crust densifying with pressure (deep_crust)."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.deep_crust import densified_mass, G_EARTH

ON = {**P.PARAMS, 'crust_density_rule': 'pressure'}


def _numeric(T, msed, p, n=20000):
    z = np.linspace((1-p['lower_crust_share'])*T, T, n+1)
    zm = .5*(z[1:]+z[:-1])
    P_ = G_EARTH*(msed+p['rho_crust']*zm)
    f = np.clip((P_-p['densify_start_gpa']*1e9)/((p['densify_neutral_gpa']-p['densify_start_gpa'])*1e9), 0, 1)
    return float(((p['rho_mantle']-p['rho_lower_crust'])*f*np.diff(z)).sum())


def test_default_is_uniform_and_gives_nothing():
    assert P.PARAMS['crust_density_rule'] == 'uniform'
    assert np.all(densified_mass(np.array([60e3]), np.array([0.]), np.array([True]), P.PARAMS) == 0)


def test_measured_values():
    assert ON['lower_crust_share'] == .324 and ON['rho_lower_crust'] == 2947.
    assert ON['densify_start_gpa'] == 1.0 and ON['densify_neutral_gpa'] == 2.0


def test_exact_integral_matches_quadrature():
    for T in (30e3, 38e3, 45e3, 55e3, 70e3, 90e3):
        for ms in (0., 2.65e6):
            m = float(densified_mass(np.array([T]), np.array([ms]), np.array([True]), ON)[0])
            assert abs(m-_numeric(T, ms, ON)) < 1e-3*max(_numeric(T, ms, ON), 1.)+1., (T, ms)


def test_shallow_crust_and_ocean_take_nothing():
    # 1 GPa at ~36 km: a 34 km column's base is shallower; oceanic nodes never densify
    m = densified_mass(np.array([34e3, 60e3]), np.zeros(2), np.array([True, False]), ON)
    assert np.all(m == 0)


def test_grows_with_thickness_and_is_held_past_neutral():
    T = np.linspace(36e3, 100e3, 65)
    m = densified_mass(T, np.zeros_like(T), np.ones_like(T, bool), ON)
    assert np.all(np.diff(m) > 0)
    # past ~72 km (2 GPa) the base is at rho_mantle; the lower crust's top deepens with the column (its share is
    # fixed: it deepens by (1 - share) per metre), so each added metre gains drho (1 - (1 - share) f(top)), f the
    # ramp's fraction at the top
    drho, s = ON['rho_mantle']-ON['rho_lower_crust'], ON['lower_crust_share']
    top = (1-s)*T[-1]
    f = np.clip((G_EARTH*ON['rho_crust']*top-1e9)/1e9, 0, 1)
    d = (m[-1]-m[-2])/(T[-1]-T[-2])
    assert abs(d-drho*(1-(1-s)*f)) < 3., (d, drho*(1-(1-s)*f))


def test_surface_sinks_by_the_load_on_a_thick_column():
    g = CubeSphere(48)
    out = {}
    for rule in ('uniform', 'pressure'):
        sim = P.PlateSim(g, seed=3, params={'crust_density_rule': rule, 'water_volume_km3': None, 'erosion': False,
                                            'threads': 1, 'lithosphere': 'none'})
        sim.state()
        w = sim.world
        cont = w['cont']
        w['thick'] = np.where(cont, 55e3, w['thick'])
        out[rule] = (sim.elevation(), cont, w['sed'].copy())
    z0, cont, sed = out['uniform']
    z1 = out['pressure'][0]
    m = densified_mass(np.full(cont.sum(), 55e3), P.sediment_column(sed, cont, P.PARAMS)[1][cont], np.ones(cont.sum(), bool), ON)
    drop = np.average((z0-z1)[cont], weights=g.area[cont])
    # the model's own flexure of that load (at n = 48 its spherical-harmonic roundtrip is not 1/rho_m: test scale)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'erosion': False, 'threads': 1})
    want = np.average(sim.flex_dry.deflection(np.where(cont, m.mean(), 0.))[cont], weights=g.area[cont])
    assert abs(drop/want-1) < .02, (drop, want)
