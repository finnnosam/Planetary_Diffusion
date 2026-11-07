"""Pass 67 (journal 89): where spontaneous subduction initiation breaks a plate that carries continent (`init_site`).
'margin' (as built): the plate's ocean component holding the old floor breaks away along the coast. 'ocean': old floor
welded to the continent does not initiate; old floor apart from it breaks off along its edge with younger floor."""
import numpy as np
import pytest

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.graphops import components

# the as-built start: continental plates whose Voronoi cells reach past the coast carry ocean
OLD_START = {'water_volume_km3': None, 'threads': 1, 'mantle': False, 'start_plates': 'voronoi',
             'start_poles': 'random', 'start_breakup': 'dispersal'}


def _mixed_plate(sim):
    c = sim.world['cont']
    best = None
    for pl in sim.plates:
        own = sim.owner == pl.id
        n_ocean = (own & ~c).sum()
        if (own & c).any() and (best is None or n_ocean > best[1]):
            best = (pl, n_ocean)
    return best[0]


def _hops_from(mask, nbr, within, n):
    """Nodes of `within` up to n graph hops from mask."""
    reach = mask.copy()
    ok = nbr >= 0
    for _ in range(n):
        reach = reach | (within & np.any(ok & reach[np.maximum(nbr, 0)], 1))
    return reach


def _setup(site, seed=3):
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=seed, params={**OLD_START, 'init_site': site})
    pl = _mixed_plate(sim)
    own = sim.owner == pl.id
    c = sim.world['cont']
    ocean = own & ~c
    coast = ocean & np.any((g.nbr >= 0) & (own & c)[np.maximum(g.nbr, 0)], 1)
    return g, sim, pl, own, ocean, coast


def test_ocean_is_the_default_and_bad_values_raise():
    assert P.PARAMS['init_site'] == 'ocean'                 # adopted pass 68 (journal 90)
    g, sim, pl, own, ocean, coast = _setup('nowhere')
    with pytest.raises(ValueError):
        sim._initiate_subduction(pl, own, coast)


def test_margin_breaks_welded_old_floor_along_the_coast():
    g, sim, pl, own, ocean, coast = _setup('margin')
    assert ocean.sum() > 50
    old = _hops_from(coast, g.nbr, ocean, 2)
    n0 = len(sim.plates)
    sim._initiate_subduction(pl, own, old)
    assert len(sim.plates) == n0+1 and sim.events[-1]['kind'] == 'subduction initiation'
    new = sim.owner == sim.events[-1]['new_plate']
    assert not (new & sim.world['cont']).any()
    assert (new & coast).sum() >= .9*coast.sum()                       # the break follows the coast


def test_ocean_leaves_welded_old_floor_and_breaks_free_old_floor_offshore():
    g, sim, pl, own, ocean, coast = _setup('ocean')
    old = _hops_from(coast, g.nbr, ocean, 2)
    n0, rng0 = len(sim.plates), sim.rng.bit_generator.state
    sim._initiate_subduction(pl, own, old)
    assert len(sim.plates) == n0 and sim.events[-1]['kind'] == 'initiation at margin skipped'
    assert sim.rng.bit_generator.state == rng0                         # no draws: the arms stay paired
    assert np.array_equal(sim.owner == pl.id, own)
    # old floor at least 4 hops offshore, not touching the continent: it breaks off, the coast stays on the plate
    near = _hops_from(coast, g.nbr, ocean, 3)
    far = ocean & ~near
    lab = components(far, g.nbr)
    big = lab == np.bincount(lab[lab >= 0]).argmax()
    assert big.sum() >= 8
    sim._initiate_subduction(pl, own, big)
    assert len(sim.plates) == n0+1 and sim.events[-1]['kind'] == 'subduction initiation'
    new = sim.owner == sim.events[-1]['new_plate']
    assert new.any() and not (new & coast).any() and not (new & sim.world['cont']).any()
    assert (sim.owner[coast] == pl.id).all()


def test_ocean_equals_margin_on_plates_without_continent():
    out = {}
    for site in ('margin', 'ocean'):
        g = CubeSphere(48)
        sim = P.PlateSim(g, seed=3, params={**OLD_START, 'init_site': site})
        c = sim.world['cont']
        pl = [p for p in sim.plates if not ((sim.owner == p.id) & c).any()][0]
        own = sim.owner == pl.id
        idx = np.flatnonzero(own)
        old = np.zeros(g.size, bool); old[idx[:idx.size//3]] = True
        sim._initiate_subduction(pl, own, old)
        out[site] = sim.owner.copy()
    assert np.array_equal(out['margin'], out['ocean'])


def test_the_rift_recorder_is_output_only_and_accounts_for_its_nodes():
    # pass 67: RiftRecorder reads PassiveRecorder's reading and draws nothing; every tracked node of a Myr is, the
    # next Myr, either still tracked (reached), become passive, or ended in one of the fates
    from planet_reviews.passive_survival import PassiveRecorder, RiftRecorder, RIFT_FATES, pool_rift
    g = CubeSphere(24)
    base = {'water_volume_km3': None, 'threads': 1}
    a = P.PlateSim(g, seed=7, params=base)
    b = P.PlateSim(g, seed=7, params=base)
    pr, rr = PassiveRecorder(g), RiftRecorder(g)
    st = b.state(); pr.update(b, st); rr.update(b, st, pr.last)
    for _ in range(25):
        a.step(); b.step()
        n_prev = 0 if rr.prev is None else rr.prev[0].shape[0]
        before = {k: dict(v) for k, v in rr.bins.items()}
        st = b.state(); pr.update(b, st); rr.update(b, st, pr.last)
        k = max(rr.bins)
        bb = rr.bins[k]; b0 = before.get(k, {f'end_{f}_km': 0. for f in RIFT_FATES})
        out = sum(bb[f'end_{f}_km']-b0.get(f'end_{f}_km', 0.) for f in RIFT_FATES)/rr.sp_km
        assert out <= n_prev+1e-6
    sa, sb = a.state(), b.state()
    for key in ('crust_thickness_m', 'elevation_m', 'plate_id'):
        assert np.array_equal(sa[key], sb[key])
    r = pool_rift([('s', rr.summary())])
    assert r is not None and 0. <= r['conversion'] <= 1.
    # particles are conserved: everything that entered either ended or is still tracked
    ent = sum(v['entered_km'] for v in rr.bins.values())
    end = sum(v[f'end_{f}_km'] for v in rr.bins.values() for f in RIFT_FATES)
    assert abs(ent-end-rr.prev[0].shape[0]*rr.sp_km) < 1e-6*max(ent, 1.)
