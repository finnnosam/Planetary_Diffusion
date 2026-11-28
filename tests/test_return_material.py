"""The material-aware return under remap_rule 'overlap' (pass 112, journal 153).

C's return thinned oceanic basement below 7 km with nothing booked (journal 150): a forearc cut on a world cell whose
overlap mixes an accreted arc with normal basement was spread over every contributor, the basement included (the
review's trace, journal 151), and the per-plate correction moved basement no process touched. The example here: a
plate's ocean at the normal 7 km with a 12 km arc strip; the world cuts the arc's excess as the forearc law does
(never more than the cell's thickness above 7 km), and erodes the continent at the same time."""
import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P

BASE = float(P.PARAMS['ocean_thick_m'])
ARC = 12e3


def _arc_sim(seed=3):
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=seed, params={'erosion': False, 'mantle': False, 'water_volume_km3': None,
                                           'remap_rule': 'overlap'})
    for _ in range(6):
        sim.state()
        sim.step()
    sim.state()
    for pl in sim.plates:
        pl.Q = fields.rotation(pl.axis, pl.speed*1e-2*1e6/g.radius)@pl.Q
    sim.t += 1.
    sim.p['hidden_steps'] = 10**9                   # nothing consumed: no collision volume in the isolated update
    k = int(np.argmax([(pl.mask & ~pl.cont).sum() for pl in sim.plates]))
    pl = sim.plates[k]
    oc = pl.mask & ~pl.cont
    pl.thick[oc] = BASE
    pl.sed[:] = 0.
    pl.rock[:] = 0.
    # the arc: a strip of the plate's ocean (a band in the plate's own frame), 12 km thick
    j = np.flatnonzero(oc)
    c = g.xyz[j].mean(0)
    c /= np.linalg.norm(c)
    s = g.xyz[j]@c
    arc = j[(s > np.quantile(s, .55)) & (s < np.quantile(s, .70))]
    pl.thick[arc] = ARC
    return sim, k, arc


def _compose(sim):
    idx, claimed = sim._claims()
    sim._world_from(idx, claimed, record=True)
    sim.new_continent = np.zeros(sim.g.size, bool)
    return idx


def _cuts(sim, k):
    """The forearc cut (at most 1.5 km, never below the cell's 7 km) on the plate's oceanic world cells and a 200 m
    subaerial cut on its continental ones; the requested volumes by crust type (world cells typed as they read)."""
    g, w = sim.g, sim.world
    ws, _, _ = sim._tent[k]
    wtype = sim._tent_type[k]
    d = np.zeros(g.size)
    oc = ws[~wtype]
    d[oc] = -np.minimum(1500., np.clip(w['thick'][oc]-BASE, 0., None))
    d[ws[wtype]] = -200.
    a = g.area
    return d, float((a[oc]*d[oc]).sum()), float((a[ws[wtype]]*d[ws[wtype]]).sum())


def _run(sim, idx, d):
    pl = [p.thick.astype(np.float64).copy() for p in sim.plates]
    ET0 = dict(sim.erosion_totals)
    sim._update_rasters(d, idx)
    diff = [p.thick.astype(np.float64)-b for p, b in zip(sim.plates, pl)]
    dET = {key: sim.erosion_totals.get(key, 0.)-ET0.get(key, 0.) for key in sim.erosion_totals}
    return pl, diff, dET


def _check(sim, k, before, diff, dET, req_oc, req_c, arc_out=()):
    g, pl, a = sim.g, sim.plates[k], sim.g.area
    oc = pl.mask & ~pl.cont
    base = oc & (before[k] == BASE)
    assert base.sum() > 100 and np.isin(np.flatnonzero(oc), np.flatnonzero(base)).mean() > .5
    assert (diff[k][base] == 0).all()                       # no basement node changes, near the arc or far from it
    assert pl.thick[oc].min() >= BASE                         # nothing below the basement
    assert dET.get('return_thick_clamped_m3', 0.) == 0
    removed_oc = float((a*diff[k])[oc].sum())
    short = dET['return_shortfall_thick_m3']                  # loss requested and not removed, booked (both types)
    assert req_oc < -1e10                                     # a real cut: > 10 km^3 requested
    # requested = removed + booked shortfall; float32 storage of a 12 km column rounds each node by ~1 mm
    tol = 2e-4*abs(req_oc)
    removed_c = float((a*diff[k])[pl.mask & pl.cont].sum())
    assert abs(req_oc+req_c-(removed_oc+removed_c)-(-short)) < tol+2e-4*abs(req_c)
    assert abs(req_c-removed_c) < 2e-4*abs(req_c)+1.          # the continental cut is conserved on its own
    for j in arc_out:                                         # nodes the plate no longer owns are left as they were
        assert diff[k][j] == 0
    return removed_oc, short


def test_forearc_cut_takes_only_arc_material_and_conserves_the_volume():
    sim, k, arc = _arc_sim()
    idx = _compose(sim)
    d, req_oc, req_c = _cuts(sim, k)
    before, diff, dET = _run(sim, idx, d)
    removed_oc, short = _check(sim, k, before, diff, dET, req_oc, req_c)
    assert short <= 1e-6*abs(req_oc)                          # every cut finds its material: nothing short
    assert abs(removed_oc-req_oc) < 2e-4*abs(req_oc)
    assert diff[k][arc].min() < -100                          # the arc itself is what thins


def test_the_old_return_would_have_thinned_the_basement():
    # the same cut spread over every contributor by area (the pre-pass-112 return, rebuilt here): basement nodes
    # beside the arc go below 7 km. Guards that the example reproduces the defect it was built for.
    sim, k, arc = _arc_sim()
    _compose(sim)
    d, _, _ = _cuts(sim, k)
    g, pl = sim.g, sim.plates[k]
    ws, cs, f = sim._tent[k]
    wt = f*g.area[ws][:, None]
    W = np.bincount(cs.ravel(), weights=wt.ravel(), minlength=g.size)
    S = np.bincount(cs.ravel(), weights=(wt*d[ws][:, None]).ravel(), minlength=g.size)
    oc = pl.mask & ~pl.cont & (pl.thick == BASE) & (W > 0)
    assert (S[oc]/W[oc] < -1.).sum() > 5


def test_ownership_change_never_moves_the_cut_onto_basement():
    # part of the arc is no longer owned by its plate (the world node over it went to another plate): the world still
    # reads it, so it takes its share of the cut like any contributor (as 'nearest' gives increments to the nodes read);
    # never onto basement, nothing short
    sim, k, arc = _arc_sim()
    idx = _compose(sim)
    d, req_oc, req_c = _cuts(sim, k)
    pl = sim.plates[k]
    wn = fields.nearest_nodes(sim.g, sim.g.xyz@pl.Q.T)
    out = arc[::3]
    other = (k+1) % len(sim.plates)
    sim.owner_slot = sim.owner_slot.copy()
    sim.owner_slot[wn[out]] = other
    before, diff, dET = _run(sim, idx, d)
    _, short = _check(sim, k, before, diff, dET, req_oc, req_c)
    assert short <= 1e-6*abs(req_oc)
    assert (diff[k][out] < 0).sum() > .5*out.size              # the unowned arc nodes lose arc material too


def _expose(sim, k, j):
    """Node j becomes one no overlap reaches (each sample that read it reads another candidate of its cell instead)."""
    ws, cs, f = sim._tent[k]
    hit = cs == j
    rows = np.flatnonzero(hit.any(1))
    for r in rows:
        alt = cs[r][cs[r] != j]
        if not alt.size:
            return False
        cs[r][cs[r] == j] = alt[0]
    return True


def test_a_newly_exposed_node_is_debited_only_against_its_own_material():
    sim, k, arc = _arc_sim()
    idx = _compose(sim)
    d, req_oc, req_c = _cuts(sim, k)
    pl = sim.plates[k]
    g = sim.g
    wn = fields.nearest_nodes(g, g.xyz@pl.Q.T)
    own = sim.owner_slot[wn] == k
    base = np.flatnonzero(pl.mask & ~pl.cont & (pl.thick == BASE) & own & ~sim.new_crust[wn])
    # a basement node beside the arc and an arc node, each exposed, each under a cut larger than its material
    nb = np.isin(g.nbr[base], arc).any(1)
    jb = next(j for j in base[nb] if _expose(sim, k, j))
    ja = next(j for j in arc if own[j] and not sim.new_crust[wn[j]] and _expose(sim, k, j))
    d[wn[jb]] = -3000.
    d[wn[ja]] = -8000.                                       # beyond the arc node's 5 km of material
    ws, _, _ = sim._tent[k]
    wtype = sim._tent_type[k]
    a = g.area
    req_oc = float((a[ws[~wtype]]*d[ws[~wtype]]).sum())
    before, diff, dET = _run(sim, idx, d)
    _check(sim, k, before, diff, dET, req_oc, req_c)
    assert diff[k][jb] == 0
    assert BASE <= pl.thick[ja] < ARC


def test_a_gain_on_a_node_the_plate_no_longer_owns_goes_to_that_node():
    # the Transfer Sweep's seed 1004 at -134 Myr (pass 112's first version): a cell's 102 m gain whose only contributor
    # was hidden, closed against the denormal gains of the two owned nodes of its type, gave an infinite factor. The
    # hidden node is what the world read, so it takes the gain; the other nodes keep their own (tiny) gains; finite
    sim, k, arc = _arc_sim()
    idx = _compose(sim)
    g, pl = sim.g, sim.plates[k]
    ws, cs, f = sim._tent[k]
    wn = fields.nearest_nodes(g, g.xyz@pl.Q.T)
    single = np.flatnonzero((cs == cs[:, :1]).all(1))         # cells read from one plate node alone
    r = single[0]
    j = cs[r, 0]
    others = np.flatnonzero(sim._tent_type[k] == sim._tent_type[k][r])
    others = others[~np.isin(cs[others], j).any(1)]
    d = np.zeros(g.size)
    d[ws[others]] = 1e-308                                    # denormal gains on every other cell of the type
    d[ws[r]] = 102.
    sim.owner_slot = sim.owner_slot.copy()
    sim.owner_slot[wn[j]] = (k+1) % len(sim.plates)           # the contributor is no longer the plate's
    before, diff, dET = _run(sim, idx, d)
    assert all(np.isfinite(p.thick).all() for p in sim.plates)
    assert abs(diff[k][j]*g.area[j]-102.*g.area[ws[r]]) < 2e-3*102.*g.area[ws[r]]   # float32 on a ~10 km column
    rest = np.ones(g.size, bool); rest[j] = False
    assert np.abs(diff[k][rest]).max() < 1e-3
    assert dET.get('return_unplaced_thick_m3', 0.) == 0
