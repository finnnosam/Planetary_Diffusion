"""C's 'lamb' frame (eightieth pass, journal 104): L by tectonothermal age measured from Hoggard et al. 2020's
seismic lithosphere on Earth's ages, Lamb's k as the density excess; the default frame unchanged."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import lithosphere as lt
from planet_tectonics import plates as P

LAMB = {**P.PARAMS, 'lith_frame': 'lamb'}


def test_default_frame_is_artemieva_and_unchanged():
    p = P.PARAMS
    assert p['lith_frame'] == 'artemieva'
    assert lt.frame(p) == (p['lith_ref_km'], p['lith_drho_kg_m3'], p['lith_max_km'])
    a = np.array([0., 50., 1000., 3600.])
    assert np.array_equal(lt.thickness_km(a, p), 93.6+.04*a)


@pytest.mark.skipif(not DATA.joinpath(*lt.HOGGARD_SL2013SV).exists() or not (DATA/'crust1'/'CNtype1-1.txt').exists(),
                    reason='Hoggard 2020 SL2013sv or CRUST1.0 not present')
def test_lamb_table_reproduces_hoggard():
    table, ref = lt.earth_lamb_table()
    assert abs(ref-lt.LAMB_REF_KM) < .01
    ages = np.array(lt.LAMB_AGE_MA)
    for t, L in table.items():
        assert abs(L-np.interp(t, ages, lt.LAMB_L_KM)) < .01, (t, L)


def test_lamb_law_monotone_invertible_and_measured_frame():
    ref, drho, top = lt.frame(LAMB)
    assert ref == lt.LAMB_REF_KM and abs(drho-35.48) < .01 and abs(top-177.7436) < 1e-6
    t = np.linspace(0., 4000., 4001)
    L = lt.thickness_km(t, LAMB)
    assert np.all(np.diff(L[:3601]) > 0) and np.all(L[3601:] == L[3600])   # strictly rising to 3.6 Ga: it inverts
    assert abs(L[0]-93.78) < 1e-9 and abs(lt.thickness_km(np.array([1750.]), LAMB)[0]-177.74) < .01
    back = lt.age_of_thickness(L[:3601], LAMB)
    assert np.allclose(back, t[:3601], atol=1e-3)
    assert lt.age_of_thickness(np.array([50.]), LAMB)[0] == 0.        # below the youngest: age 0
    dz = lt.elevation_term(np.array([0., 3200., np.nan]), LAMB)
    assert abs(dz[0]-lt.LAMB_K_M_PER_KM*(ref-93.78)) < .01 and dz[2] == 0.   # Lamb's k: 10.753 m per km of L
    assert abs(dz[1]+lt.LAMB_K_M_PER_KM*(177.74-ref)) < .1
    s = lt.stretched_age(np.array([1750.]), np.array([.5]), LAMB)       # 177.7 -> 88.9 km: below the youngest
    assert s[0] == 0.


def test_lamb_thickening_bounded_at_the_table_top():
    therm = np.array([100.-1750., 100.-33., np.nan])
    th, m = lt.thickened_therm(therm, np.full(3, 40e3), np.array([20e3, 4e3, 5e3]), 100., LAMB)
    L = lt.thickness_km(100.-th[:2], LAMB)
    assert abs(L[0]-lt.frame(LAMB)[2]) < 1e-6                           # held at the cratonic top
    assert abs(L[1]-lt.thickness_km(np.array([33.]), LAMB)[0]*1.1) < 1e-6 and np.isnan(th[2]) and not m[2]


def test_lamb_frame_runs():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'lith_frame': 'lamb', 'threads': 1})
    sim.state()
    for _ in range(12):
        sim.step(); st = sim.state()
    c = st['continental']; la = st['lithosphere_age_myr']
    assert not np.isnan(la[c]).any()
    L = lt.thickness_km(la[c], LAMB)
    assert L.min() >= 93.78-1e-9 and L.max() <= lt.frame(LAMB)[2]+1e-6
