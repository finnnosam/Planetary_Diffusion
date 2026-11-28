"""Loose sediment is entrained without bedrock's incision threshold (sediment_entrainment 'cover'; journal 110)."""
import numpy as np
import pytest

import planet_tectonics.plates as P
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl


def _setup():
    g = CubeSphere(24)
    x, y, zc = g.xyz.T
    # a low-gradient continent (a plateau falling gently to the sea) and a bumpy deep floor
    z = np.where(x > .3, 600.*(x-.3)+400.+30.*np.cos(9*y), -4000.+2000.*np.cos(6*y)*np.cos(5*zc))
    p = {**gl.DEFAULTS, 'open_edges': False, 'theta': 200.}
    return g, z, p, z < -1500., np.full(g.size, gl.DEFAULTS['K'])


def _step(g, z, p, seed, K, cover=None, dt=2.5e5):
    return gl.surface_step(g, z, seed, None, np.ones(g.size, bool), K, p, dt, g.spacing_min(), g.diffusion_bound(),
                           keep_deep=True, cover=cover)


def test_no_cover_is_the_bedrock_law_exactly():
    g, z, p, seed, K = _setup()
    r0 = _step(g, z, p, seed, K)
    r1 = _step(g, z, p, seed, K, cover=np.zeros(g.size))
    assert np.array_equal(r0['z'], r1['z']) and np.array_equal(r0['deep_m3'], r1['deep_m3'])


def test_cover_erodes_below_the_threshold_and_never_below_its_base():
    g, z, p, seed, K = _setup()
    land = z > 0
    r0 = _step(g, z, p, seed, K)
    cov = np.where(land, 50., 0.)
    r = _step(g, z, p, seed, K, cover=cov)
    e0 = (np.maximum(z-r0['z'], 0)*g.area)[land].sum()
    e = (np.maximum(z-r['z'], 0)*g.area)[land].sum()
    assert e > 1.5*e0                                        # the low-gradient land loses its loose cover
    # where the bedrock law alone does not cut, the surface falls no deeper than the cover's base (hillslope
    # diffusion aside: compare the fluvial kernel, which surface_step applies before diffusion)
    zf = gl.priority_flood(z, gl.connected_ocean(z, seed, g.nbr), np.ones(g.size, bool), g.nbr, 1e-3)
    rec, dist, stack = gl.route(zf, gl.connected_ocean(z, seed, g.nbr), g.nbr, g.dist)
    area = gl.accumulate(rec, stack, g.area)
    ocean = gl.connected_ocean(z, seed, g.nbr)
    a0 = gl.erode_and_deposit(z, zf, ocean, rec, dist, stack, area, K, p['m'], p['G'], p['Sc'], 2.5e5, g.area,
                              p['lake_tolerance'], 6, p['theta'])[0]
    a1 = gl.erode_and_deposit_cover(z, zf, ocean, rec, dist, stack, area, K, p['m'], p['G'], p['Sc'], 2.5e5, g.area,
                                    p['lake_tolerance'], 6, p['theta'], cov)[0]
    bed = a0 >= z-1e-9                                       # the bedrock law does not erode here
    assert (a1[bed] >= z[bed]-cov[bed]-1e-6).all()
    assert (z-a1)[land & bed].max() > 1.                     # and somewhere it does take cover


def test_cover_keeps_the_volume_closed():
    g, z, p, seed, K = _setup()
    r = _step(g, z, p, seed, K, cover=np.where(z > 0, 80., 20.))
    eroded = (r['eroded']*g.area).sum()
    assert abs(((r['z']-z)*g.area).sum()+r['deep_m3'].sum()) < 1e-6*eroded


def test_plate_layer_cover_exports_more_and_keeps_the_ledger():
    g = CubeSphere(48)
    out = {}
    for rule in ('bedrock', 'cover'):
        # (pass 112, journal 153) pins remap_rule 'nearest' / rift_area 'world': with C (overlap, plates) as the default this
        # trajectory redeposits 0.853 per unit eroded against
        # 0.808 (journal 150 sec. 3); the subject is the cover rule's direction on the state it was set on
        sim = P.PlateSim(g, seed=5, params={'sediment_entrainment': rule, 'threads': 1, 'remap_rule': 'nearest',
                                            'rift_area': 'world',   # burial needs 'cover'
                                            'sediment_rock': 'burial' if rule == 'cover' else 'none'})
        sim.state()
        v0 = sim.continental_volume()
        for _ in range(12):
            sim.step()
            sim.state()
        t = sim.erosion_totals
        out[rule] = t
        assert np.isfinite(sim.world['thick']).all() and (sim.world['sed'] >= -1e-6).all()
        assert abs(sim.continental_volume()-v0)/v0 < .05
    # loose deposits leave the land: less is redeposited on continents per unit eroded
    keep = {r: out[r]['continental_deposited_m3']/-out[r]['continental_eroded_m3'] for r in out}
    assert keep['cover'] < keep['bedrock']


def test_default_is_cover_and_unknown_rule_raises():
    assert P.PARAMS['sediment_entrainment'] == 'cover'
    sim = P.PlateSim(CubeSphere(16), seed=1, params={'sediment_entrainment': 'nope', 'threads': 1})
    sim.state()
    with pytest.raises(ValueError):
        sim.step()
