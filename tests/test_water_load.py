"""Pass 94 (journal 124): water_load_rule. 'parts' (as built) water-loads a thinned crust's own Airy depth wherever it
lies below sea level, also under a sediment fill booked as a dry load; 'column' decides once per column, by its
air-loaded local level (1-D Airy), as elevation_check's gate reads it."""
import numpy as np
import pytest

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.isostasy import sediment_column

BASE = {'water_volume_km3': None, 'erosion': False, 'threads': 1, 'lithosphere': 'none',
        'crust_density_rule': 'uniform'}


def _surface(rule, solid_sed_m, thick_m=31e3):
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={**BASE, 'water_load_rule': rule})
    sim.state()
    w = sim.world
    cont = w['cont']
    w['thick'] = np.where(cont, thick_m, w['thick'])
    w['sed'] = np.where(cont, solid_sed_m, w['sed'])
    sim.sea_level = 0.
    return sim, sim.elevation(), cont


def test_default_is_column():
    assert P.PARAMS['water_load_rule'] == 'column'      # adopted pass 95 (journal 126)
    with pytest.raises(ValueError):
        _surface('nonsense', 0.)


def _levels(sim, cont, solid):
    p = sim.p
    rc, rm, rw = p['rho_crust'], p['rho_mantle'], p['rho_water']
    sh, sm_dry, _ = sediment_column(np.where(cont, solid, 0.), cont, p)
    own = p['cont_ref_elev_m']+(31e3-p['cont_ref_thick_m'])*(rm-rc)/rm
    return own, sh, sm_dry, rm, rw


def test_dry_fill_is_not_water_loaded_under_column():
    solid = 2500.                                     # fills the thinned crust's depth above sea level
    s0, z0, cont = _surface('parts', solid)
    s1, z1, _ = _surface('column', solid)
    own, sh, sm_dry, rm, rw = _levels(s1, cont, solid)
    assert own < 0
    assert np.all(own+sh[cont]-sm_dry[cont]/rm > 0)                  # dry on the column's own level
    assert np.all(own*rm/(rm-rw)+sh[cont] > 0)                       # and on the as-built criterion
    # the same dry load either way: the rules differ by the crust's water term alone, on every continental node
    d = (z0-z1)[cont]
    assert np.allclose(d, own*rw/(rm-rw), atol=1e-6), (d.min(), d.max(), own*rw/(rm-rw))
    assert np.allclose(z0[~cont], z1[~cont], atol=1e-6)


def test_wet_column_is_the_same_under_both():
    solid = 50.                                        # a thin cover: the column stays under the sea
    s0, z0, cont = _surface('parts', solid)
    s1, z1, _ = _surface('column', solid)
    own, sh, sm_dry, rm, rw = _levels(s1, cont, solid)
    assert np.all(own+sh[cont]-sm_dry[cont]/rm < 0) and np.all(own*rm/(rm-rw)+sh[cont] < 0)
    assert np.allclose(z0, z1, atol=1e-6)
