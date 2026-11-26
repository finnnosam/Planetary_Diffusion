"""remap_rule (pass 108, journal 145): how values move between the world grid and the rotated plate rasters.

'nearest' (as built) copies the plate node nearest each world node's position in the owner's frame and returns a step's
increments to the plate nodes read. 'tent' interpolates from the nearest plate node and its eight neighbours of the same
crust type (exact at the nodes, a partition of unity) and returns the increments by the same weights, conserved per
plate. The floor's second trace (journal 145) found the composition's per-step jump in the plate rasters' departure from
the world, not in the resampling of the world's own surface."""
import numpy as np
import pytest

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P


def _moved_sim(remap, seed=3):
    """A small run whose plates have moved, then one more rotation and composition, ready for an isolated
    `_update_rasters` call with chosen increments (as tests/test_return_transfer.py)."""
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=seed, params={'erosion': False, 'mantle': False, 'water_volume_km3': None,
                                           'remap_rule': remap})
    for _ in range(6):
        sim.state()
        sim.step()
    sim.state()
    for pl in sim.plates:
        pl.Q = fields.rotation(pl.axis, pl.speed*1e-2*1e6/g.radius)@pl.Q
    sim.t += 1.
    idx, claimed = sim._claims()
    sim._world_from(idx, claimed, record=True)
    sim.new_continent = np.zeros(g.size, bool)
    sim.p['hidden_steps'] = 10**9                   # nothing consumed: no collision volume in the isolated update
    return sim, idx


def _owned_live(sim, k):
    g, pl = sim.g, sim.plates[k]
    wn = fields.nearest_nodes(g, g.xyz@pl.Q.T)
    J = np.flatnonzero(pl.mask & (sim.owner_slot[wn] == k) & ~sim.new_crust[wn])
    return J, wn[J]


def _biggest(sim):
    return int(np.argmax([pl.mask.sum() for pl in sim.plates]))


def _apply(sim, idx, dthick, dsed=None):
    before = [pl.thick.astype(np.float64).copy() for pl in sim.plates]
    sb = [pl.sed.astype(np.float64).copy() for pl in sim.plates]
    sim._update_rasters(dthick, idx, dsed, None)
    return ([pl.thick.astype(np.float64)-b for pl, b in zip(sim.plates, before)],
            [pl.sed.astype(np.float64)-b for pl, b in zip(sim.plates, sb)])


def test_default_is_nearest_and_bad_settings_raise():
    assert P.PARAMS['remap_rule'] == 'nearest'
    sim, idx = _moved_sim('nearest')
    sim.p['remap_rule'] = 'bilinear'
    with pytest.raises(ValueError):
        sim._world_from(*sim._claims())
    sim, idx = _moved_sim('tent')
    sim.p['return_transfer'] = 'conservative'
    with pytest.raises(ValueError):
        sim._update_rasters(np.zeros(sim.g.size), idx)


def test_weights_exact_at_nodes_and_a_partition_of_unity():
    sim, idx = _moved_sim('tent')
    g = sim.g
    k = _biggest(sim)
    pl = sim.plates[k]
    ws, cs, wt = sim._tent[k]
    assert ws.size > 100
    assert np.allclose(wt.sum(1), 1.) and (wt >= 0).all()
    assert pl.mask[cs].all()
    assert (pl.cont[cs][wt > 0] == np.repeat(pl.cont[cs[:, 0]][:, None], 9, 1)[wt > 0]).all()   # one crust type
    # a world position exactly on a plate node reads that node alone
    Q = pl.Q
    pl.Q = np.eye(3)
    j = np.flatnonzero(pl.mask)[:500]
    c2, w2 = sim._tent_weights(pl, j, j)
    pl.Q = Q
    assert np.allclose(w2[:, 0], 1.) and np.allclose(w2[:, 1:], 0.)
    assert (c2[:, 0] == j).all()


def test_overlap_weights_exact_at_nodes_and_sum_near_one_per_plate_node():
    # 'overlap': a world cell exactly over a plate cell reads it alone; on a moved plate every plate node's overlap
    # weights (in area) sum to about its own area, so the return keeps a uniform increment uniform without rescaling
    sim, idx = _moved_sim('overlap')
    g = sim.g
    k = _biggest(sim)
    pl = sim.plates[k]
    ws, cs, wt = sim._tent[k]
    assert np.allclose(wt.sum(1), 1.)
    W = np.bincount(cs.ravel(), weights=(wt*g.area[ws][:, None]).ravel(), minlength=g.size)
    J, _ = _owned_live(sim, k)
    interior = J[(W[J] > 0)]
    r = W[interior]/g.area[interior]
    assert np.median(np.abs(r-1.)) < .15
    Q = pl.Q
    pl.Q = np.eye(3)
    j = np.flatnonzero(pl.mask)[:500]
    c2, w2 = sim._overlap_weights(pl, j, j)
    pl.Q = Q
    assert (c2 == j[:, None]).all()


@pytest.mark.parametrize('rule', ['tent', 'overlap'])
def test_interpolation_follows_a_smooth_field_better_than_nearest(rule):
    # a smooth field in the plate's frame on its continental nodes: the tent read is closer to the field at each world
    # node's position than the nearest node's copy
    sim, idx = _moved_sim(rule)
    g = sim.g
    k = max(range(len(sim.plates)), key=lambda i: (sim.plates[i].mask & sim.plates[i].cont).sum())
    pl = sim.plates[k]
    f = lambda xyz: 35e3+3e3*np.sin(6*xyz[:, 0])*np.cos(5*xyz[:, 1])
    pl.thick[:] = f(g.xyz).astype(np.float32)
    idx, claimed = sim._claims()
    sim._world_from(idx, claimed, record=False)
    ws, cs, wt = sim._tent[k]
    c = sim.world['cont'][ws]
    truth = f(g.xyz[ws]@pl.Q)[c]
    tent = sim.world['thick'][ws][c]
    near = pl.thick[idx[k, ws]].astype(np.float64)[c]
    e_t, e_n = np.sqrt(np.mean((tent-truth)**2)), np.sqrt(np.mean((near-truth)**2))
    assert e_t < .6*e_n


@pytest.mark.parametrize('rule', ['tent', 'overlap'])
def test_uniform_increment_leaves_every_live_owned_node_the_same(rule):
    # a uniform 100 m increment on every world node: every live owned node of the plate gets the same increment within
    # its crust type (none is left behind), and each type's gain is the world's over the plate's nodes of that type
    sim, idx = _moved_sim(rule)
    g = sim.g
    k = _biggest(sim)
    J, _ = _owned_live(sim, k)
    ws = np.flatnonzero((sim.owner_slot == k) & ~sim.new_crust)
    cw = sim.world['cont'].copy()
    pc = sim.plates[k].cont.copy()
    d, _ = _apply(sim, idx, np.full(g.size, 100.))
    for typ in (True, False):                        # each crust type uniform, its volume the world's of that type
        Jt, wt = J[pc[J] == typ], ws[cw[ws] == typ]
        if not Jt.size:
            continue
        assert np.ptp(d[k][Jt]) < .05
        assert abs(np.median(d[k][Jt])-100.) < 10.
        assert abs(float((g.area[Jt]*d[k][Jt]).sum())/float((g.area[wt]*100.).sum())-1.) < 1e-4


@pytest.mark.parametrize('rule', ['tent', 'overlap'])
def test_smooth_increment_field_stays_smooth(rule):
    def rough(g, J, d):
        inJ = np.zeros(g.size, bool); inJ[J] = True
        nb = g.nbr[J]
        ok = (nb >= 0) & inJ[np.where(nb >= 0, nb, 0)]
        m = ok.any(1)
        nbv = np.where(ok, d[np.where(nb >= 0, nb, 0)], 0.)
        return float(np.mean(np.abs(d[J][m]-nbv[m].sum(1)/ok[m].sum(1))))
    sim, idx = _moved_sim(rule)
    g = sim.g
    k = _biggest(sim)
    J, wJ = _owned_live(sim, k)
    field = 100.+50.*(g.xyz@np.array([.6, .48, .64]))
    d, _ = _apply(sim, idx, field)
    ref = np.zeros(g.size); ref[J] = field[wJ]
    assert rough(g, J, d[k]) <= 1.05*rough(g, J, ref)+.05
    assert np.corrcoef(d[k][J], field[wJ])[0, 1] > .999


@pytest.mark.parametrize('rule', ['tent', 'overlap'])
def test_sediment_losses_and_gains_are_conserved_per_plate(rule):
    sim, idx = _moved_sim(rule)
    g = sim.g
    k = _biggest(sim)
    J, _ = _owned_live(sim, k)
    pl = sim.plates[k]
    pl.sed[:] = 5000.
    sim.world['sed'][:] = 5000.
    side = g.xyz@np.array([0., 0., 1.])
    ds = np.where(side > 0, 20.*side, -10.*np.abs(side))
    ws = np.flatnonzero((sim.owner_slot == k) & ~sim.new_crust)
    _, d = _apply(sim, idx, np.zeros(g.size), ds)
    a = g.area
    gain_w = float((a[ws]*np.maximum(ds[ws], 0)).sum()); loss_w = float((a[ws]*np.minimum(ds[ws], 0)).sum())
    gain_p = float((a[J]*np.maximum(d[k][J], 0)).sum()); loss_p = float((a[J]*np.minimum(d[k][J], 0)).sum())
    assert abs(gain_p/gain_w-1.) < 1e-3
    assert abs(loss_p/loss_w-1.) < 1e-3


@pytest.mark.parametrize('rule', ['tent', 'overlap'])
def test_tent_run_keeps_both_ledgers(rule):
    # the whole step with erosion and continental sediment: the crust and continental-sediment volumes change only by
    # their booked causes, as under 'nearest' (test_plates)
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'remap_rule': rule})
    v0, s0 = sim.continental_volume(), sim.continental_sediment_volume()
    for _ in range(15):
        sim.state()
        sim.step()
    v, s1 = sim.continental_volume(), sim.continental_sediment_volume()
    assert abs((v-v0)-sum(sim.ledger.values())) < 2e-3*v0
    assert abs((s1-s0)-sum(sim.sed_ledger.values())) < 1e-6*max(s1, 1.)+1e-3
    assert sim.erosion_totals['return_rescaled_thick_m3'] > 0
    st = sim.state()
    assert np.isfinite(st['elevation_m']).all() and (st['sediment_m'] >= 0).all()
