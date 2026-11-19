"""start_sediment (eighty-ninth pass, journal 115): the start's continents carry Earth's stable cover, the start's floor
its own pelagic rain."""
from pathlib import Path

import numpy as np
import pytest

from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import plates as P
from planet_tectonics.isostasy import sediment_column
from planet_tectonics.start import _nscore


def test_start_sediment_all_is_the_default():
    # adopted 17 Nov (ninetieth pass, journal 116); 'none' leaves the start bare
    assert P.PARAMS['start_sediment'] == 'all'
    sim = P.PlateSim(CubeSphere(48), seed=3, params={'erosion': False, 'start_sediment': 'none'})
    assert not sim.world['sed'].any()


def test_start_sediment_leaves_the_rest_of_the_start_unchanged():
    g = CubeSphere(48)
    a = P.PlateSim(g, seed=3, params={'threads': 1, 'start_sediment': 'none'})
    b = P.PlateSim(g, seed=3, params={'threads': 1, 'start_sediment': 'all'})
    for k in ('cont', 'thick', 'born', 'orog', 'therm'):
        assert np.array_equal(a.world[k], b.world[k], equal_nan=True), k
    assert np.array_equal(a.owner, b.owner)
    for x, y in zip(a.plates, b.plates):                 # the plates' stream is untouched (own stream, seed 9)
        assert np.array_equal(x.axis, y.axis) and x.speed == y.speed
    assert b.world['sed'][b.world['cont']].sum() > 0


def test_start_cover_takes_earths_distribution_and_correlation():
    g = CubeSphere(96)
    sim = P.PlateSim(g, seed=5, params={'start_sediment': 'continent'})
    w = sim.world; c = w['cont']; a = g.area
    h = sediment_column(w['sed'], c, sim.p)[0][c]
    o = np.argsort(h); cw = np.cumsum(a[c][o]); cw /= cw[-1]
    for q in (.5, .75, .9):                              # bulk cover at Earth's quantiles (interpolated table)
        want = np.interp(q, P.EARTH_STABLE_CRUST_Q, P.EARTH_STABLE_SED_BULK_M)
        assert abs(h[o][np.searchsorted(cw, q)]-want) < .1*want+50
    mean = float((h*a[c]).sum()/a[c].sum())
    assert 1150 < mean < 1500                            # Earth's stable cover: 1,331 m bulk
    assert not w['sed'][~c].any()                        # 'continent' leaves the floor bare
    inner = c & (gl.distance_to(~c, g.nbr, g.dist, 5000e3) > 200e3)   # off the thinned rim
    zs, zt = _nscore(w['sed'][inner].astype(float), a[inner]), _nscore(w['thick'][inner].astype(float), a[inner])
    r = float(np.sum(a[inner]*zs*zt)/np.sqrt(np.sum(a[inner]*zs**2)*np.sum(a[inner]*zt**2)))
    assert -.45 < r < -.15                               # Earth's -0.31: thin crust under thick cover
    # coherent basins, not node noise: neighbours' scores correlate as Earth's cover does at < 100 km (0.84)
    z = np.zeros(g.size); z[c] = _nscore(w['sed'][c].astype(float), a[c])
    nb = g.nbr; ok = (nb >= 0) & c[:, None] & c[np.clip(nb, 0, None)]
    i, j = np.nonzero(ok)
    assert np.corrcoef(z[i], z[nb[i, j]])[0, 1] > .7


def test_start_floor_carries_its_own_pelagic_rain():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'start_sediment': 'all'})
    w = sim.world; o = ~w['cont']
    age = sim.t-w['born'][o].astype(float)
    assert np.allclose(w['sed'][o], sim.p['pelagic_m_myr']*age, rtol=1e-5, atol=1e-3)


def test_start_sediment_keeps_the_ledgers():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'start_sediment': 'all', 'threads': 1})
    v0, s0 = sim.continental_volume(), sim.continental_sediment_volume()
    assert s0 > 0
    for _ in range(12):
        sim.step()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
    s1 = sim.continental_sediment_volume()
    assert abs((s1-s0)-sum(sim.sed_ledger.values())) < 1e-6*max(s1, 1.)+1e-3
    st = sim.state()
    assert np.isfinite(st['elevation_m']).all() and (st['sediment_m'] >= 0).all()


def test_earth_stable_sediment_matches_crust1():
    """EARTH_STABLE_SED_BULK_M and EARTH_STABLE_SED_RHO are CRUST1.0's on EARTH_STABLE_CRUST_KM's domain."""
    root = Path(__file__).resolve().parents[1]
    if not (root/'data'/'crust1'/'crsthk.xyz').exists():
        pytest.skip('CRUST1.0 not present')
    from planet_reviews.crust_review import wpct
    from planet_reviews.interior_probe import GROUPS
    from planet_earth.bird import earth_boundaries
    from planet_earth.crust1 import C1, crust1_cell
    g = CubeSphere(217); a = g.area; cell = crust1_cell(g)
    T = np.loadtxt(C1/'crsthk.xyz')[:, 2][cell]*1e3; S = np.loadtxt(C1/'sedthk.xyz')[:, 2][cell]*1e3
    ty = np.array([l.split() for l in (C1/'CNtype1-1.txt').read_text().splitlines() if l.strip()]).ravel()[cell]
    dbnd = gl.distance_to(earth_boundaries(g), g.nbr, g.dist, 5000e3)/1e3
    m = (T >= 20e3) & np.isin(ty, GROUPS['Precambrian/platform']) & (dbnd > 500)
    got = [wpct(S[m], a[m], 100*q) for q in P.EARTH_STABLE_CRUST_Q]
    assert np.allclose(got, P.EARTH_STABLE_SED_BULK_M, atol=1.)
    zs, zt = _nscore(S[m], a[m]), _nscore((T-S)[m], a[m])
    r = np.sum(a[m]*zs*zt)/np.sqrt(np.sum(a[m]*zs**2)*np.sum(a[m]*zt**2))
    assert abs(r-P.EARTH_STABLE_SED_RHO) < 1e-3
