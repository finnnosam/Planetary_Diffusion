"""return_transfer (pass 107, journal 144): how a step's world increments return to the plate rasters.

'read' (as built) gives a world node's increment to the plate node it read, so plate nodes no world node reads miss the
step (the trace of passes 105-106, journal 143). 'conservative' gives every live plate node the plate still owns the
increment of the world node over it, with each plate's gains and losses rescaled to the world volume of the nodes it
owns."""
import numpy as np
import pytest

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P


def _moved_sim(rule, seed=3, compose='nearest'):
    """A small run whose plates have moved (so the forward and backward maps are not inverses), then one more
    rotation and composition, ready for an isolated `_update_rasters` call with chosen increments."""
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=seed, params={'erosion': False, 'mantle': False, 'water_volume_km3': None, 'remap_rule': 'nearest',
                                           'return_transfer': rule, 'compose_map': compose})
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
    """Plate k's live nodes whose world location it owns (and that are not this step's new crust), and the world
    node over each."""
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


def test_default_is_read_and_unknown_rule_raises():
    assert P.PARAMS['return_transfer'] == 'read'
    sim, idx = _moved_sim('read')
    sim.p['return_transfer'] = 'nearest'
    with pytest.raises(ValueError):
        sim._update_rasters(np.zeros(sim.g.size), idx)


def test_uniform_increment_leaves_a_rigid_plate_uniform():
    # a uniform 100 m increment on every world node: under 'conservative' every live owned node of the plate gets the
    # same increment, and the plate's volume gain is the world's over the nodes it owns; under 'read' some nodes no
    # world node reads get nothing (the defect the trace found)
    out = {}
    for rule in ('read', 'conservative', 'matched'):
        sim, idx = _moved_sim(rule)
        g = sim.g
        k = _biggest(sim)
        J, _ = _owned_live(sim, k)
        ws = np.flatnonzero((sim.owner_slot == k) & ~sim.new_crust)
        d, _ = _apply(sim, idx, np.full(g.size, 100.))
        out[rule] = (d[k][J], float((g.area[ws]*100.).sum()), float((g.area[J]*d[k][J]).sum()))
    for rule in ('conservative', 'matched'):
        dc, world, plate = out[rule]
        assert np.ptp(dc) < .05                      # float32 rounding on ~35 km of crust
        assert abs(np.median(dc)-100.) < 10.         # the area rescaling is a few per cent
        assert abs(plate/world-1.) < 1e-4 or rule == 'matched'   # matched also feeds read nodes the plate hides
    dr = out['read'][0]
    assert dr.min() < 1.                             # 'read' leaves unread nodes behind


def test_smooth_increment_field_stays_smooth():
    # a smooth increment field (a large-scale gradient across the globe): under 'conservative' the plate's increments
    # are no rougher than the field itself sampled on the plate; 'read' adds the unread nodes' zeros and the
    # mean-of-two readers
    def rough(sim, k, J, d):
        g = sim.g
        inJ = np.zeros(g.size, bool); inJ[J] = True
        nb = g.nbr[J]
        ok = (nb >= 0) & inJ[np.where(nb >= 0, nb, 0)]
        m = ok.any(1)
        nbv = np.where(ok, d[np.where(nb >= 0, nb, 0)], 0.)
        return float(np.mean(np.abs(d[J][m]-nbv[m].sum(1)/ok[m].sum(1))))
    res = {}
    for rule in ('read', 'conservative'):
        sim, idx = _moved_sim(rule)
        g = sim.g
        k = _biggest(sim)
        J, wJ = _owned_live(sim, k)
        field = 100.+50.*(g.xyz@np.array([.6, .48, .64]))
        d, _ = _apply(sim, idx, field)
        ref = np.zeros(g.size); ref[J] = field[wJ]
        res[rule] = (rough(sim, k, J, d[k]), rough(sim, k, J, ref), np.corrcoef(d[k][J], field[wJ])[0, 1])
    rc, rref, corr = res['conservative']
    assert rc <= 1.05*rref+.05
    assert corr > .999
    assert res['read'][0] > 2*rref


def test_sediment_losses_and_gains_are_conserved_per_plate():
    # erosion on half the globe, deposition on the other: each part's plate volume equals the world volume of the
    # nodes the plate owns (where no sediment runs out: the losses are kept below the cover)
    sim, idx = _moved_sim('conservative')
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


@pytest.mark.parametrize('rule,compose', [('conservative', 'nearest'), ('matched', 'nearest'), ('state', 'nearest'),
                                          ('read', 'matched')])
def test_conservative_run_keeps_both_ledgers(rule, compose):
    # the whole step with erosion and continental sediment: the crust and continental-sediment volumes change only
    # by their booked causes, as under 'read' (test_plates)
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'return_transfer': rule, 'compose_map': compose,
                                        'remap_rule': 'nearest'})
    v0, s0 = sim.continental_volume(), sim.continental_sediment_volume()
    for _ in range(15):
        sim.state()
        sim.step()
    v, s1 = sim.continental_volume(), sim.continental_sediment_volume()
    assert abs((v-v0)-sum(sim.ledger.values())) < 2e-3*v0
    assert abs((s1-s0)-sum(sim.sed_ledger.values())) < 1e-6*max(s1, 1.)+1e-3
    if rule in ('conservative', 'matched'):
        assert sim.erosion_totals['return_rescaled_thick_m3'] > 0
    st = sim.state()
    assert np.isfinite(st['elevation_m']).all() and (st['sediment_m'] >= 0).all()


def test_matched_composition_reads_fewer_shared_nodes():
    # compose_map 'matched': ownership as 'nearest', then world nodes that share a plate node with a nearer world node
    # move to a free live node within two rings of it: fewer shared reads, never a node outside the plate's mask, and
    # never further than two rings and a half from the world node
    res = {}
    for compose in ('nearest', 'matched'):
        sim, idx = _moved_sim('read', compose=compose)
        g = sim.g
        shared, far = 0, 0.
        for k, pl in enumerate(sim.plates):
            ws = np.flatnonzero((sim.owner_slot == k) & ~sim.new_crust)
            js = idx[k, ws]
            assert pl.mask[js].all()
            shared += int((np.bincount(js, minlength=g.size) > 1).sum())
            far = max(far, float(np.sqrt(((g.xyz[ws]@pl.Q-g.xyz[js])**2).sum(1)).max()*g.radius))
        res[compose] = (shared, far)
    assert res['matched'][0] < res['nearest'][0]
    assert res['matched'][1] < 3.5*sim.h_min*1.6
    sim.p['compose_map'] = 'bijective'
    with pytest.raises(ValueError):
        sim._world_from(*sim._claims())
