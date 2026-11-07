"""Pass 66 (journal 88): the start's breakup (`start_breakup` 'whole': the supercontinent's plates as one plate) and
the passive margins' survival record (planet_reviews.passive_survival, output only)."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_reviews.passive_survival import PassiveRecorder, FATES, pool, passive_edge

BASE = {'water_volume_km3': None, 'threads': 1}


def _cont_plates(sim):
    c = sim.world['cont']
    return sorted({int(i) for i in sim.owner[c]})


def test_whole_is_the_default_and_puts_the_supercontinent_on_one_plate():
    # 'whole' (with start_plates 'coast', start_poles 'forces') is the default since the sixty-seventh pass (journal 89);
    # compared here on the as-built Voronoi start with random poles
    assert P.PARAMS['start_breakup'] == 'whole'
    g = CubeSphere(24)
    old = {**BASE, 'start_plates': 'voronoi', 'start_poles': 'random'}
    a = P.PlateSim(g, seed=5, params={**old, 'start_breakup': 'dispersal'})
    b = P.PlateSim(g, seed=5, params={**old, 'start_breakup': 'whole'})
    assert np.array_equal(a.world['cont'], b.world['cont'])          # the same supercontinent
    assert np.array_equal(a.world['thick'], b.world['thick'])
    assert len(_cont_plates(a)) >= 5 and len(_cont_plates(b)) == 1     # as built, oceanic plates reach past the coast
    assert len(b.plates) == len(a.plates)-4
    # the oceanic plates draw the same numbers: same poles and speeds, in the same order
    oa, ob = a.plates[5:], b.plates[1:]                               # the oceanic seeds' plates, in creation order
    assert [round(p.speed, 12) for p in oa] == [round(p.speed, 12) for p in ob]
    assert all(np.allclose(p.axis, q.axis) for p, q in zip(oa, ob))
    pg = [pl for pl in b.plates if pl.cont.any()][0]
    assert 2. <= pg.speed <= 4.
    for _ in range(3):
        b.step()


def test_whole_with_the_coast_start_ends_the_continental_plate_at_the_coast():
    g = CubeSphere(24)
    b = P.PlateSim(g, seed=5, params={**BASE, 'start_breakup': 'whole', 'start_plates': 'coast'})
    pid = _cont_plates(b)
    assert len(pid) == 1
    assert not ((b.owner == pid[0]) & ~b.world['cont']).any()


def test_the_recorder_is_output_only_and_accounts_for_every_passive_node():
    g = CubeSphere(24)
    a = P.PlateSim(g, seed=7, params=BASE)
    b = P.PlateSim(g, seed=7, params=BASE)
    rec = PassiveRecorder(g)
    rec.update(b, b.state())
    n_prev = []
    for _ in range(4):
        a.step(); b.step()
        n_prev.append(rec.prev[0].shape[0])
        st = b.state()
        before = {k: dict(v) if isinstance(v, dict) else v for k, v in rec.bins.items()}
        rec.update(b, st)
    sa, sb = a.state(), b.state()
    for k in ('crust_thickness_m', 'elevation_m', 'plate_id', 'sediment_m'):
        assert np.array_equal(sa[k], sb[k])
    # every node passive last Myr either survived or ended in one of the fates
    s = rec.summary()
    ended = sum(bn[f'end_{f}_km'] for bn in s['bins'] for f in FATES)/s['node_km']
    assert ended >= 0 and np.isfinite(ended)
    r = pool([('seed7', s)])
    assert r is not None and r['passive_km'] >= 0 and abs(sum(r['fates'].values())-(1 if ended > 0 else 0)) < 1e-9


def test_a_still_planet_keeps_its_passive_margins():
    """Plates that do not move: the record carries every passive node onto itself, so nothing ends."""
    g = CubeSphere(24)
    sim = P.PlateSim(g, seed=7, params=BASE)
    for pl in sim.plates:
        pl.speed = 0.
    st = sim.state()
    st['boundary_type'] = np.zeros(g.size, np.int8)                  # no boundaries: every edge node is passive
    rec = PassiveRecorder(g)
    rec.update(sim, st)
    rec.update(sim, st)
    pas = passive_edge(st, g, sim.p)[0]
    b = rec.summary()['bins'][-1]
    assert pas.sum() > 0
    assert sum(b[f'end_{f}_km'] for f in FATES) == 0 and b['births_rifted_km'] == b['births_converted_km'] == 0
