"""rift_area (pass 110, journal 148): the new area a pure-shear rift gap opens. 'world' (as built) strains the zone by the
gap's world cell; 'plates' by the area the owning plate's raster switches on over the gap, so the zone gives exactly the
volume the new plate nodes take (the nearest-node maps between the world and a rotated raster are not one to one)."""
import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from test_plates import _rift_line


def test_rift_area_default_is_plates():
    assert P.PARAMS['rift_area'] == 'plates'             # since pass 112 (journal 153)


def test_kernel_opened_area_conserves_and_defaults_to_the_cell():
    g = CubeSphere(217)
    gap, zone, _ = _rift_line(g)
    thick = np.where(zone, 38e3, 0.)
    rows = np.flatnonzero(gap)
    base = P.rift_pure_shear(rows, g.area, thick, zone, g.nbr, g.dist, 100e3, 10e3)
    same = P.rift_pure_shear(rows, g.area, thick, zone, g.nbr, g.dist, 100e3, 10e3, g.area.copy())
    for x, y in zip(base, same):                       # the cell's own area: the same arithmetic, bit for bit
        assert np.array_equal(x, y)
    rng = np.random.default_rng(1)
    opened = g.area*rng.uniform(.3, 1.8, g.size)
    opened[rows[::5]] = 0.                             # some gaps open nothing on the plates
    Tbar, Tg, dT = P.rift_pure_shear(rows, g.area, thick, zone, g.nbr, g.dist, 100e3, 10e3, opened)
    assert (Tg[::5] == 0).all()
    # the zone loses exactly what the gaps receive on the area they open
    assert abs((dT*g.area).sum()+(Tg*opened[rows]).sum()) < 1e-9*(Tg*opened[rows]).sum()
    assert (dT[~zone] == 0).all() and (dT <= 0).all()
    # less opened area in total, less thinning
    assert -(dT*g.area).sum() < -(base[2]*g.area).sum()*opened[rows].sum()/g.area[rows].sum()*1.05


def test_opened_area_is_what_the_update_switches_on():
    """`_opened_area` at rift time is, per world node, the area of the plate nodes `_update_rasters` switches on over new
    crust (holes the update fills afterwards, copied from a neighbour and booked as 'holes', can only add to it)."""
    g = CubeSphere(48)
    p = {'erosion': False, 'rift_rule': 'pure_shear', 'rift_width_km': 600., 'mantle': False, 'rift_area': 'plates',
         'remap_rule': 'overlap', 'subduction_H_cont_m': 0., 'arc_H_ocean_m': 0., 'subduction_erosion_H_m': 0.}
    sim = P.PlateSim(g, seed=3, params=p)
    o_up, o_rps = sim._update_rasters, sim._rift_pure_shear
    seen = {'steps': 0, 'match': 0, 'gaps': 0}

    def rps(*a):
        seen['A'] = sim._opened_area()                  # at rift time
        return o_rps(*a)

    def up(dthick, idx, dsed=None, drock=None):
        A = seen['A']
        before = [pl.mask.copy() for pl in sim.plates]
        wns = [fields.nearest_nodes(g, g.xyz@pl.Q.T) for pl in sim.plates]
        nc = sim.new_crust.copy()
        out = o_up(dthick, idx, dsed, drock)
        got = np.zeros(g.size)
        for pl, m0, wn in zip(sim.plates, before, wns):
            on = pl.mask & ~m0 & nc[wn]
            np.add.at(got, wn[on], g.area[on])
        assert (got >= A*(1-1e-12)).all()
        k = A > 0
        seen['gaps'] += int(k.sum())
        seen['match'] += int(np.isclose(got[k], A[k], rtol=1e-12).sum())
        seen['steps'] += 1
        return out

    sim._rift_pure_shear, sim._update_rasters = rps, up
    for _ in range(25):
        seen['A'] = np.zeros(g.size)
        sim.step()
    assert seen['steps'] == 25 and seen['gaps'] > 0
    assert seen['match'] >= .95*seen['gaps']


def test_plates_rule_closes_the_rift_on_the_plates():
    """One step with a rift: the plates' continental volume change from new rift nodes equals the zone's thinning on the
    plates, under 'overlap' (exact return) and 'plates' (exact opening), within the gaps' own later increments."""
    g = CubeSphere(48)
    p = {'erosion': False, 'rift_rule': 'pure_shear', 'rift_width_km': 600., 'mantle': False, 'rift_area': 'plates',
         'remap_rule': 'overlap', 'crust_flow_rule': 'fitted', 'crust_kappa_normal_m2_yr': 0., 'crust_kappa_m2_yr': 0.,
         'subduction_H_cont_m': 0., 'arc_H_ocean_m': 0., 'subduction_erosion_H_m': 0., 'margin_shortening': False,
         'backarc_extension': False}
    sim = P.PlateSim(g, seed=3, params=p)
    o_rps = sim._rift_pure_shear
    rec = {}

    def rps(*a):
        c0 = sim.world['cont'].copy()
        out = o_rps(*a)
        r = np.flatnonzero(sim.world['cont'] & ~c0 & sim.new_crust)
        A = sim._opened_area()
        rec['zone'] = -float((out*g.area).sum())
        rec['gap'] = float((A[r]*sim.world['thick'][r]).sum())
        return out

    sim._rift_pure_shear = rps
    n = 0
    for _ in range(40):
        rec.clear()
        sim.step()
        if rec.get('zone', 0.) > 0:
            n += 1
            assert abs(rec['gap']-rec['zone']) <= 1e-9*rec['zone']
    assert n > 0
