"""start_cover (ninety-sixth pass, journal 129): the start's cover lies where its column has room for it
('accommodation'), with Earth's correlations of the stable cover with the crust and the bare level."""
from pathlib import Path

import numpy as np
import pytest

from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import plates as P
from planet_tectonics.isostasy import sediment_column
from planet_tectonics.lithosphere import elevation_term
from planet_tectonics.start import _nscore


def _corr(a, x, y):
    return float(np.sum(a*x*y)/np.sqrt(np.sum(a*x*x)*np.sum(a*y*y)))


def _scores(sim):
    """Normal scores of the cover, the crust and the bare level on the continent off the thinned rim."""
    g, w, p = sim.g, sim.world, sim.p
    c = w['cont']
    inner = c & (gl.distance_to(~c, g.nbr, g.dist, 5000e3) > 200e3)
    T = w['thick'].astype(float)
    bare = (p['cont_ref_elev_m']+(T-p['cont_ref_thick_m'])*(1.-p['rho_crust']/p['rho_mantle'])
            + sim.mantle.topography(sim.t)+elevation_term(sim.t-w['therm'].astype(float), p))
    a = g.area[inner]
    return a, _nscore(w['sed'][inner].astype(float), a), _nscore(T[inner], a), _nscore(bare[inner], a)


def test_default_is_accommodation_and_unknown_rules_raise():
    # 'accommodation' adopted 20 Nov (ninety-seventh pass, journal 130)
    assert P.PARAMS['start_cover'] == 'accommodation'
    with pytest.raises(ValueError):
        P.PlateSim(CubeSphere(48), seed=3, params={'start_cover': 'basins'})
    bare = P.PlateSim(CubeSphere(48), seed=3, params={'start_cover': 'accommodation', 'start_sediment': 'none'})
    assert not bare.world['sed'][bare.world['cont']].any()       # a bare start: start_cover is moot (journal 130)


def test_accommodation_moves_only_the_cover():
    g = CubeSphere(48)
    a = P.PlateSim(g, seed=3, params={'threads': 1, 'start_cover': 'copula'})
    b = P.PlateSim(g, seed=3, params={'threads': 1, 'start_cover': 'accommodation'})
    for k in ('cont', 'thick', 'born', 'orog', 'therm'):
        assert np.array_equal(a.world[k], b.world[k], equal_nan=True), k
    assert np.array_equal(a.owner, b.owner)
    for x, y in zip(a.plates, b.plates):
        assert np.array_equal(x.axis, y.axis) and x.speed == y.speed
    o = ~a.world['cont']
    assert np.array_equal(a.world['sed'][o], b.world['sed'][o])          # the floor's pelagic rain unchanged
    c = a.world['cont']
    assert not np.array_equal(a.world['sed'][c], b.world['sed'][c])
    ha = sediment_column(a.world['sed'], c, a.p)[0][c]; hb = sediment_column(b.world['sed'], c, b.p)[0][c]
    area = g.area[c]
    for q in (.5, .75, .9):                                               # the same distribution: ranks moved only
        qa = ha[np.argsort(ha)][np.searchsorted(np.cumsum(area[np.argsort(ha)])/area.sum(), q)]
        qb = hb[np.argsort(hb)][np.searchsorted(np.cumsum(area[np.argsort(hb)])/area.sum(), q)]
        assert abs(qa-qb) < .1*qa+50
    assert abs((ha*area).sum()-(hb*area).sum()) < .05*(ha*area).sum()


def test_accommodation_takes_earths_correlations():
    g = CubeSphere(96)
    res = {}
    for rule in ('copula', 'accommodation'):
        sim = P.PlateSim(g, seed=5, params={'start_cover': rule})
        a, zs, zt, zb = _scores(sim)
        r_st, r_sb, r_tb = _corr(a, zs, zt), _corr(a, zs, zb), _corr(a, zt, zb)
        res[rule] = (r_st, r_sb, (r_sb-r_st*r_tb)/np.sqrt((1-r_st**2)*(1-r_tb**2)))
    r_st, r_sb, part = res['accommodation']
    assert -.5 < r_st < -.15                             # Earth's -0.333 (the rim's thinning off: inner nodes)
    assert -.9 < r_sb < -.6                              # Earth's -0.787
    assert part < -.5 and part < res['copula'][2]-.4     # Earth's partial -0.76; the copula's is near zero


def test_accommodation_keeps_the_ledgers():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'start_cover': 'accommodation', 'threads': 1})
    v0, s0 = sim.continental_volume(), sim.continental_sediment_volume()
    for _ in range(12):
        sim.step()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
    s1 = sim.continental_sediment_volume()
    assert abs((s1-s0)-sum(sim.sed_ledger.values())) < 1e-6*max(s1, 1.)+1e-3
    st = sim.state()
    assert np.isfinite(st['elevation_m']).all() and (st['sediment_m'] >= 0).all()


def test_earth_stable_bare_correlations_match_crust1():
    """EARTH_STABLE_SED_BARE is CRUST1.0 / ETOPO's on the stable domain without the ice sheets."""
    root = Path(__file__).resolve().parents[1]
    if not (root/'data'/'crust1'/'crsthk.xyz').exists():
        pytest.skip('CRUST1.0 not present')
    from planet_grid.paths import OUTPUTS
    if not (OUTPUTS/'planet-reference'/'earth-n217.npz').exists():
        pytest.skip('earth-n217.npz not present')
    from planet_reviews.interior_probe import GROUPS
    from planet_earth.bird import earth_boundaries
    from planet_earth.crust1 import C1, crust1_cell
    from planet_earth.reference import ice_mask
    g = CubeSphere(217); a = g.area; cell = crust1_cell(g); p = P.PARAMS
    T = np.loadtxt(C1/'crsthk.xyz')[:, 2][cell]*1e3; S = np.loadtxt(C1/'sedthk.xyz')[:, 2][cell]*1e3
    ty = np.array([l.split() for l in (C1/'CNtype1-1.txt').read_text().splitlines() if l.strip()]).ravel()[cell]
    dbnd = gl.distance_to(earth_boundaries(g), g.nbr, g.dist, 5000e3)/1e3
    z = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
    m = (T >= 20e3) & np.isin(ty, GROUPS['Precambrian/platform']) & (dbnd > 500) & ~ice_mask(g)
    b = z-S*(1-p['rho_sediment']/p['rho_mantle'])
    zs, zt, zb = _nscore(S[m], a[m]), _nscore((T-S)[m], a[m]), _nscore(b[m], a[m])
    assert abs(_corr(a[m], zs, zt)-P.EARTH_STABLE_SED_BARE[0]) < 1e-3
    assert abs(_corr(a[m], zs, zb)-P.EARTH_STABLE_SED_BARE[1]) < 1e-3
