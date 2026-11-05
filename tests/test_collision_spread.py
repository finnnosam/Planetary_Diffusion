"""Pass 64 (journal 86): the collision spread (`collision_spread` 'earth') and R, the eclogite foundering
(`eclogite_rule` 'foundering'). Both default off; the default path is unchanged."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.crust_kernels import COLLISION_KERNEL, COLLISION_KERNEL_KM, collision_kernel


def test_new_rules_are_off_by_default():
    assert P.PARAMS['collision_spread'] == 'node'
    assert P.PARAMS['eclogite_rule'] == 'none'
    assert P.PARAMS['eclogite_thick_m'] == 65000. and P.PARAMS['eclogite_tau_myr'] == 10.


def test_collision_kernel_is_earths_table_and_zero_beyond_it():
    assert collision_kernel(0.) == 1. and collision_kernel(150e3) == COLLISION_KERNEL[1]
    assert collision_kernel(COLLISION_KERNEL_KM[-1]*1e3+1.) == 0.
    assert np.all(np.diff(COLLISION_KERNEL[:5]) < 0)          # densest at the suture


def test_spread_conserves_volume_and_stays_on_the_receiving_plate():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'collision_spread': 'earth'})
    sim.step()
    w = sim.world
    cont = np.flatnonzero(w['cont'])
    i = cont[np.argmax(np.bincount(sim.owner[cont]).take(sim.owner[cont]))]   # a node on the most continental plate
    dV = np.zeros(g.size); dV[i] = 1e15
    out = sim._spread_collision(dV)
    assert abs(out.sum()-1e15) < 1e-6*1e15
    got = out > 0
    assert got.sum() > 10 and out[i] < 1e15
    assert np.all(sim.owner[got] == sim.owner[i]) and np.all(w['cont'][got] | (np.arange(g.size) == i)[got])
    # density per unit area falls with distance, as Earth's kernel
    from planet_landscape import graph_landscape as gl
    src = np.zeros(g.size, bool); src[i] = True
    d = gl.distance_to(src, g.nbr, g.dist, 4e6)
    near, far = got & (d < 200e3), got & (d > 600e3)
    if far.any():
        assert (out[near]/g.area[near]).mean() > (out[far]/g.area[far]).mean()


def test_spread_and_foundering_runs_keep_the_ledger():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'collision_spread': 'earth', 'eclogite_rule': 'foundering',
                                        'eclogite_thick_m': 40000., 'crust_kappa_m2_yr': 0.})
    v0 = sim.continental_volume()
    for _ in range(15):
        sim.step()
    assert sim.erosion_totals.get('eclogite_foundered_m3', 0.) > 0
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0


def test_foundering_removes_the_excess_with_its_time_constant():
    from planet_tectonics.crust_kernels import eclogite_loss
    th = np.array([40e3, 65e3, 80e3, 120e3, 80e3])
    cont = np.array([1, 1, 1, 1, 0], bool)
    d = eclogite_loss(th, cont, 65e3, 10., 1.)
    assert d[0] == 0 and d[1] == 0 and d[4] == 0                  # thin, at the threshold, oceanic: untouched
    assert np.allclose(th[2:4]+d[2:4], 65e3+(th[2:4]-65e3)*np.exp(-.1))
    after = th.copy()
    for _ in range(100):
        after = after+eclogite_loss(after, cont, 65e3, 10., 1.)
    assert np.allclose(after[2:4], 65e3+(th[2:4]-65e3)*np.exp(-10.))  # e-folding of the excess in 10 Myr
