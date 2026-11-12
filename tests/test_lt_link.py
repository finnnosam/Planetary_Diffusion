"""Pass 81 (journal 105): the lithosphere joined to its crust. start_lithosphere 'crust' (the start's ages coupled to
its own crust within each distance bin by Earth's measured rho) and belt_lithosphere 'shorten' (shortened arc-belt
crust thickens its lithosphere by pure shear instead of the arc's reset); both off by default."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import lithosphere as lt
from planet_tectonics import plates as P


def test_defaults_off():
    assert P.PARAMS['start_lithosphere'] == 'distance' and P.PARAMS['belt_lithosphere'] == 'reset'


@pytest.mark.skipif(not DATA.joinpath(*lt.HOGGARD_SL2013SV).exists() or not (DATA/'crust1'/'CNtype1-1.txt').exists(),
                    reason='Hoggard 2020 SL2013sv or CRUST1.0 not present')
def test_start_rho_reproduces_earth():
    rho = lt.earth_start_rho()
    assert np.allclose(rho, lt.EARTH_START_RHO, atol=5e-4)
    assert all(0 <= r < .5 for r in rho)                                # weak, positive in every distance bin


def test_copula_join_rho0_is_the_distance_join_and_rho_couples_age_to_crust():
    rng = np.random.default_rng(0)
    n = 4000
    d = rng.uniform(0., 3000e3, n); T = rng.normal(38e3, 4e3, n); a = rng.uniform(.5, 1.5, n)
    v = np.repeat(rng.uniform(size=8), n//8)
    zero = lt.start_ages_crust(d, T, a, v, rho=np.zeros(len(lt.EARTH_START_RHO)))
    assert np.allclose(zero, lt.start_ages(d, v), atol=1e-3)            # rho 0: the craton rank alone
    full = lt.start_ages_crust(d, T, a, v, rho=np.full(len(lt.EARTH_START_RHO), .9))
    base = lt.start_ages_crust(d, T, a, v)
    tab = np.array(lt.EARTH_LITH_AGE_TABLE)
    for ages in (full, base):
        assert ages.min() >= tab.min() and ages.max() <= tab.max()
    c = lambda x: np.corrcoef(T, lt.thickness_km(x, {'lith_frame': 'lamb'}))[0, 1]
    assert c(full) > c(base) > c(zero)


def test_start_crust_join_raises_the_start_L_T_link():
    g = CubeSphere(96)
    out = {}
    for rule in ('distance', 'crust'):
        sim = P.PlateSim(g, seed=1012, params={'threads': 1, 'start_lithosphere': rule})
        st = sim.state()
        T = st['crust_thickness_m']; la = st['lithosphere_age_myr']; m = st['continental'] & np.isfinite(la)
        L = lt.thickness_km(la[m], sim.p)
        out[rule] = (np.corrcoef(T[m], L)[0, 1], L.mean(), st['elevation_m'])
    assert out['crust'][0] > out['distance'][0]
    assert abs(out['crust'][1]-out['distance'][1]) < 5.                # the same ages, re-ranked within each bin
    assert np.all(np.isfinite(out['crust'][2]))
    with pytest.raises(ValueError):
        P.PlateSim(CubeSphere(24), seed=1, params={'start_lithosphere': 'nope'})


def test_shortened_belt_nodes_are_not_reset():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'belt_lithosphere': 'shorten', 'threads': 1})
    sim.state()
    cont = np.flatnonzero(sim.world['cont'])
    sim.belt = np.zeros(g.size, bool); sim.belt[cont[:10]] = True
    sim.belt_gain = np.zeros(g.size); sim.belt_gain[cont[:4]] = 500.
    r = sim._thermal_resets(np.zeros(g.size, bool))
    assert r.sum() == 6 and not r[cont[:4]].any() and r[cont[4:10]].all()
    sim.p = {**sim.p, 'belt_lithosphere': 'reset'}
    assert sim._thermal_resets(np.zeros(g.size, bool)).sum() == 10


def test_shorten_runs_bounded_and_thickens_belt_lithosphere():
    g = CubeSphere(48)
    out = {}
    for rule in ('reset', 'shorten'):
        sim = P.PlateSim(g, seed=5, params={'belt_lithosphere': rule, 'threads': 1})
        sim.state()
        gained = np.zeros(g.size, bool)
        for _ in range(40):
            sim.step()
            st = sim.state()
            gained |= getattr(sim, 'belt_gain', np.zeros(g.size)) > 0
        out[rule] = (st, gained)
    st, gained = out['shorten']
    assert gained.any()                                                  # margin shortening reached a belt
    c = st['continental'] & np.isfinite(st['lithosphere_age_myr'])
    L = lt.thickness_km(st['lithosphere_age_myr'][c], P.PARAMS)
    assert L.max() <= lt.frame(P.PARAMS)[2]+.04*40+1e-6 and np.all(np.isfinite(st['elevation_m']))
    def L_young_thick(s):
        m = s['continental'] & (s['crust_thickness_m'] > 42e3) & (s['lithosphere_age_myr'] < 250)
        return lt.thickness_km(s['lithosphere_age_myr'][m], P.PARAMS).mean() if m.any() else np.nan
    a, b = L_young_thick(out['reset'][0]), L_young_thick(st)
    assert np.isnan(a) or np.isnan(b) or b >= a-1.
    with pytest.raises(ValueError):
        sim = P.PlateSim(g, seed=5, params={'belt_lithosphere': 'nope'}); sim.step()
