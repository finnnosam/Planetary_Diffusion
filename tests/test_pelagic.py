"""Pass 74 (journal 97; default since pass 75, journal 99): pelagic_m_myr, sediment rained on the ocean floor (Earth's far-ocean rate, 1.1 m/Myr of grain).
It lands on oceanic nodes only, its ledger term is what it adds, and the ocean-floor sediment budget closes with it."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P


def _sim(rate, seed=3):
    g = CubeSphere(48)
    # plumbing size; no erosion, so the ocean floor's only source is the pelagic rain; sea level held at the datum
    return P.PlateSim(g, seed=seed, params={'erosion': False, 'water_volume_km3': None, 'pelagic_m_myr': rate,
                                            'threads': 1})


def test_default_is_earths_rate():
    # adopted in pass 75 (journal 99, Pelagic Sweep)
    assert P.PARAMS['pelagic_m_myr'] == 1.1


def test_pelagic_lands_on_the_ocean_floor_and_its_budget_closes():
    sim = _sim(1.1)
    sim.state()
    g = sim.g
    sed0 = sim.state()['sediment_m'].astype(np.float64)
    cont0 = sim.state()['continental'].astype(bool)
    sim.step()
    st = sim.state()
    added = sim.erosion_totals['pelagic_m3']
    ocean_area = g.area[~cont0].sum()
    # one step of 1 Myr: 1.1 m on (nearly) the ocean's area; the ocean area moves a little within the step
    assert abs(added/(1.1*ocean_area)-1) < .05
    d = st['sediment_m'].astype(np.float64)-sed0
    # nothing lands on continental crust that stayed continental
    keep = cont0 & st['continental'].astype(bool)
    assert np.abs(d[keep]).max() < 1e-6
    for _ in range(9):
        sim.step()
        sim.state()
    e = sim.erosion_totals
    stock = (g.area*sim.state()['sediment_m'].astype(np.float64)).sum()-(g.area*sed0).sum()
    booked = e['pelagic_m3']-e.get('sediment_subducted_m3', 0.)-e.get('sediment_accreted_m3', 0.)-e.get(
        'collision_sediment_to_crust_m3', 0.)-e.get('basin_to_continent_m3', 0.)
    # plate-raster resampling moves sediment between nodes without conserving it exactly (STATE: ~17% of increments on
    # continents at n = 217): the books close to a few percent of what was added
    assert abs(stock-booked) < .1*e['pelagic_m3'], (stock, booked, e['pelagic_m3'])


def test_off_adds_nothing():
    sim = _sim(0.)
    sim.state()
    sim.step()
    assert sim.erosion_totals.get('pelagic_m3', 0.) == 0.
