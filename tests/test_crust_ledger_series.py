"""The stats series carries the plate layer's crust ledger (pass 95, journal 126): crust_budget's residual held merges',
pieces', holes' and the raster update's volume, which the series did not record."""
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from scripts import simulate_plates as S


def test_series_carries_the_crust_ledger_and_it_closes():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'threads': 1})
    a = S.summary(sim, sim.state())
    for _ in range(6):
        sim.step()
        st = sim.state()
    b = S.summary(sim, st)
    assert set(sim.ledger) == set(b['crust_ledger_km3'])
    dv = b['continental_volume_km3']-a['continental_volume_km3']
    dl = sum(b['crust_ledger_km3'].values())-sum(a['crust_ledger_km3'].values())
    assert abs(dv-dl) < 2e-3*a['continental_volume_km3']
