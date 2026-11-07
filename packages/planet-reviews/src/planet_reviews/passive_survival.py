"""Passive margins' survival through a run (sixty-sixth pass, journal 88): what ends them, and when; and since the
sixty-seventh pass (journal 89) the young rifted margin's record (`RiftRecorder`): how much of it becomes passive edge.

  python -m planet_reviews.passive_survival <sweep or combination or seed folders...> [--json out.json]

Recorded in the run (`PassiveRecorder`, called every Myr by scripts.simulate_plates after the step; output only: it
reads the state and draws nothing, so a run is the same with or without it). Passive edge as the scorecard reads it
at 46 km: label-free continent (crust + sediment >= 20 km) beside a node that is not, more than 300 km from any plate
boundary (planet_reviews.margin_width.BOUNDARY_KM). Each Myr every passive-edge node is carried with its plate over
the next step (the rotation the step will apply: the plate's pole and speed as they stand after the step's events)
to the node it lands on. It survives if that node or a neighbour is continental edge with no boundary within KEEP_KM
(200 km: hysteresis, so that boundary nodes wiggling by a node around 300 km do not end and rebirth margins; one node
of slack for the nearest-node raster); its birth time goes to the passive nodes there (the oldest wins). Otherwise it
ended, by what is at the landing node after the step (boundary types within 300 km):
  subduction / collision / ridge / transform: still continental edge, now within 300 km of a boundary of that type
                                               (the nearest one): the margin turned active, or a rift reopened it;
  closed:   still continental, no longer beside an ocean node (its ocean closed or filled), no boundary within 300 km;
  closed_collision: the same with a collision boundary within 300 km;
  lost:     not continental at the landing node and its neighbours (crust consumed, eroded under 20 km, or lost);
  other:    continental, not edge, with a non-collision boundary within 300 km.
Passive-edge nodes no survivor reached are births, by the oldest ocean floor beside them: < 25 Myr 'rifted' (the
ridge has just moved 300 km away), else 'converted' (a boundary nearby stopped or moved off: an active margin turned
passive). Nodes passive at the start carry birth = the start (their age is unknown: 'start' cohort).

Earth: Bradley (2008, Earth-Sci. Rev. 91, 1-26): 76 ancient passive margins have a mean lifespan of 181 Myr (130 Myr
for those of Permian to Neogene age); today's have a mean age of 104 Myr (maximum 180). In a steady state, mean
lifespan = passive-edge length / the length that ends per Myr; Earth's 173k km at 46 km (margin_width) with 130-181
Myr means ~960-1,330 km of passive edge ending per Myr. The planet's lifespan is the tracked stock (passive edge and
the nodes the hysteresis keeps) / the length ending per Myr.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_tectonics.params import RIDGE, SUBDUCTION, COLLISION, TRANSFORM

BOUNDARY_KM = 300.
KEEP_KM = 200.          # hysteresis: a passive node lives on while the nearest boundary stays beyond this
BIN_MYR = 10.
FATES = ('subduction', 'collision', 'ridge', 'transform', 'closed', 'closed_collision', 'other', 'lost')
TYPE_NAME = {RIDGE: 'ridge', SUBDUCTION: 'subduction', COLLISION: 'collision', TRANSFORM: 'transform'}
AGE_EDGES = (0., 5., 25., 50., 100., 150., 200., np.inf)       # lifespan at death, Myr
RIFTED_MYR = 25.
EARTH_LIFESPAN_MYR = (130., 181.)                                # Bradley 2008: Permian-Neogene / all ancient margins
EARTH_PASSIVE_KM = 173.1e3                                       # margin_width at 46 km (journal 85)


def passive_edge(st, g, params):
    """(passive edge, label-free continent, edge, boundary distance per type, nearest boundary distance) on the world
    nodes of a state."""
    from planet_landscape import graph_landscape as gl
    from planet_reviews.crust_review import column_thickness_m, label_free_continent
    nb = g.nbr; ok = nb >= 0; nbc = np.where(ok, nb, 0)
    cont = label_free_continent(column_thickness_m(st, params))
    edge = cont & np.any(np.where(ok, ~cont[nbc], False), 1)
    bt = np.asarray(st['boundary_type'])
    cut = 1.2*BOUNDARY_KM*1e3
    dist = {k: gl.distance_to(bt == k, nb, g.dist, cut) for k in TYPE_NAME}
    dmin = np.minimum.reduce(list(dist.values()))
    return edge & (dmin > BOUNDARY_KM*1e3), cont, edge, dist, dmin


class PassiveRecorder:
    """Every Myr: passive-edge nodes carried with their plates, their survival, ends and births (see the module)."""

    def __init__(self, g):
        self.g = g
        self.sp_km = float(np.sqrt(g.area.mean()))/1e3
        self.bins = {}
        self.prev = None                                   # (positions, birth, start cohort) of last Myr's passive nodes

    def _bin(self, t):
        k = float(np.floor(t/BIN_MYR)*BIN_MYR)
        if k not in self.bins:
            self.bins[k] = {'t0_myr': k, 'steps': 0, 'passive_km': 0., 'tracked_km': 0., 'births_rifted_km': 0., 'births_converted_km': 0.,
                            **{f'end_{f}_km': 0. for f in FATES},
                            **{f'end_{f}_lasting_km': 0. for f in FATES},
                            'end_by_age_km': [0.]*(len(AGE_EDGES)-1), 'end_age_sum_km_myr': 0., 'end_start_cohort_km': 0.}
        return self.bins[k]

    def update(self, sim, st):
        from planet_grid import fields
        g = self.g
        nb = g.nbr; ok = nb >= 0; nbc = np.where(ok, nb, 0)
        pas, cont, edge, dist, dmin = passive_edge(st, g, sim.p)
        self.last = (pas, cont, edge, dist, dmin)          # RiftRecorder reads the same Myr's reading
        b = self._bin(sim.t)
        b['steps'] += 1
        b['passive_km'] += float(pas.sum())*self.sp_km
        birth = np.full(g.size, np.inf)
        start = np.zeros(g.size, bool)
        reached = np.zeros(g.size, bool)
        if self.prev is not None:
            x, b0, s0 = self.prev
            j = fields.nearest_nodes(g, x)
            hood = np.concatenate([j[:, None], np.where(ok[j], nb[j], j[:, None])], 1)    # landing node + ring
            keep_ok = edge & (dmin > KEEP_KM*1e3)
            alive = keep_ok[hood].any(1)
            # survivors hand their birth to the passive nodes of their landing ring (the oldest wins)
            hs, hb = hood[alive].ravel(), np.repeat(b0[alive], hood.shape[1])
            hst = np.repeat(s0[alive], hood.shape[1])
            keep = keep_ok[hs]
            np.minimum.at(birth, hs[keep], hb[keep])
            reached[hs[keep]] = True
            np.logical_or.at(start, hs[keep], hst[keep])
            dead = ~alive
            jd = j[dead]; hd = hood[dead]
            near = dmin[jd] <= BOUNDARY_KM*1e3
            types = np.array(list(TYPE_NAME))
            nearest = types[np.argmin(np.stack([dist[k][jd] for k in types]), 0)] if jd.size else np.zeros(0, int)
            is_edge = edge[hd].any(1); is_cont = cont[hd].any(1)
            fate = np.full(jd.size, 'other', dtype=object)
            for k, name in TYPE_NAME.items():
                fate[is_cont & is_edge & near & (nearest == k)] = name
            fate[is_cont & ~is_edge & ~near] = 'closed'
            fate[is_cont & ~is_edge & near & (nearest == COLLISION)] = 'closed_collision'
            fate[~is_cont] = 'lost'
            age = sim.t-b0[dead]
            for f in FATES:
                b[f'end_{f}_km'] += float((fate == f).sum())*self.sp_km
            fresh = ~s0[dead]
            lasting = ~fresh | (age >= AGE_EDGES[1])
            for f in FATES:
                b[f'end_{f}_lasting_km'] += float(((fate == f) & lasting).sum())*self.sp_km
            h, _ = np.histogram(age[fresh], bins=AGE_EDGES)
            b['end_by_age_km'] = [u+float(v)*self.sp_km for u, v in zip(b['end_by_age_km'], h)]
            b['end_age_sum_km_myr'] += float(age[fresh].sum())*self.sp_km
            b['end_start_cohort_km'] += float((~fresh).sum())*self.sp_km
        birth = np.where(pas | reached, birth, np.inf)
        new = pas & ~reached
        if self.prev is not None and new.any():
            oage = np.where(cont, -np.inf, np.nan_to_num(np.asarray(st['ocean_age_myr'], np.float64), nan=-np.inf))
            nage = np.where(ok, oage[nbc], -np.inf).max(1)[new]
            b['births_rifted_km'] += float((nage < RIFTED_MYR).sum())*self.sp_km
            b['births_converted_km'] += float((nage >= RIFTED_MYR).sum())*self.sp_km
        birth[new] = sim.t
        start[new] = self.prev is None
        # carry this Myr's passive nodes with their plates over the next step
        I = np.flatnonzero(pas | (reached & np.isfinite(birth)))   # passive, and kept by the hysteresis
        b['tracked_km'] += float(I.size)*self.sp_km
        x = g.xyz[I].copy()
        rot = {pl.id: fields.rotation(pl.axis, pl.speed*1e-2*sim.p['dt_myr']*1e6/g.radius) for pl in sim.plates}
        own = sim.owner[I]
        for pid in np.unique(own):
            if pid in rot:
                m = own == pid
                x[m] = x[m]@rot[pid].T
        self.prev = (x, birth[I], start[I])

    def summary(self):
        return {'bin_myr': BIN_MYR, 'node_km': self.sp_km, 'age_edges_myr': [float(a) for a in AGE_EDGES],
                'bins': [self.bins[k] for k in sorted(self.bins)]}


RIFT_YOUNG_MYR = 10.     # the ocean floor beside a rifted margin's edge when it enters the rift record
RIFT_STALL_MYR = 100.    # a rifted margin still beside a ridge or transform this long ends as 'stalled'
RIFT_FATES = ('passive', 'subduction', 'collision', 'closed', 'lost', 'other', 'stalled')


class RiftRecorder:
    """Young rifted margin until it becomes passive edge or ends (sixty-seventh pass, journal 89; output only).

    Read with the continental flag, which carries the thinned rifted rim (the label-free reading, crust + sediment
    >= 20 km, sits inland of it and missed most rifts: seed 1008, 181 km/Myr entered against ~2,600 of rifted
    production, a quarter 'lost'). Enters: flagged continental edge with a ridge the nearest boundary within 300 km
    and ocean floor younger than RIFT_YOUNG_MYR on it or its ring, not already tracked. Each Myr every tracked node is carried with its plate over the next
    step (as PassiveRecorder) and read at its landing node and ring: 'passive' when the label-free passive edge
    (the scorecard's) lies within two rings (the margin is born: the time since entry is recorded); kept while the
    ring holds flagged continental edge whose nearest
    boundary within 300 km is a ridge or a transform (the young ocean still opening); else it ends: 'subduction' /
    'collision' (edge, the nearest boundary within 300 km of that type: the young ocean closed or a trench reached
    it), 'closed' (continent, no longer edge), 'lost' (no continent), 'other'; and 'stalled' when still opening after
    RIFT_STALL_MYR. Earth reads no such record; its steady state (rifted production 1,553 km/Myr by
    planet_reviews.rift_production's instrument on Mueller et al. 2019, against ~960-1,330 km/Myr of passive edge
    ending, Bradley 2008) puts ~0.6-0.85 of its rifted margin into passive edge."""

    def __init__(self, g):
        self.g = g
        self.sp_km = float(np.sqrt(g.area.mean()))/1e3
        self.bins = {}
        self.prev = None                                   # (positions, entry time)

    def _bin(self, t):
        k = float(np.floor(t/BIN_MYR)*BIN_MYR)
        if k not in self.bins:
            self.bins[k] = {'t0_myr': k, 'steps': 0, 'tracked_km': 0., 'entered_km': 0.,
                            **{f'end_{f}_km': 0. for f in RIFT_FATES}, 'passive_age_sum_km_myr': 0.}
        return self.bins[k]

    def update(self, sim, st, reading):
        """Each tracked piece of margin is one particle (one node's length) carried with the plate of the node it lands
        on, so entered = ended + the change of the tracked stock; the ring is read for its state, and the opening
        nodes of a kept particle's ring do not enter again."""
        from planet_grid import fields
        g = self.g
        nb = g.nbr; ok = nb >= 0; nbc = np.where(ok, nb, 0)
        pas, _, _, dist, dmin = reading
        cont = np.asarray(st['continental']).astype(bool)          # the flag: it carries the thinned rim
        edge = cont & np.any(np.where(ok, ~cont[nbc], False), 1)
        pas2 = pas | np.any(np.where(ok, pas[nbc], False), 1)      # passive edge within one ring of the node
        b = self._bin(sim.t)
        b['steps'] += 1
        near = dmin <= BOUNDARY_KM*1e3
        types = np.array(list(TYPE_NAME))
        nearest = types[np.argmin(np.stack([dist[k] for k in types]), 0)]
        opening = edge & near & ((nearest == RIDGE) | (nearest == TRANSFORM))
        covered = np.zeros(g.size, bool)
        xs, es, js = [np.zeros((0, 3))], [np.zeros(0)], [np.zeros(0, np.int64)]
        if self.prev is not None:
            x, e0 = self.prev
            j = fields.nearest_nodes(g, x)
            hood = np.concatenate([j[:, None], np.where(ok[j], nb[j], j[:, None])], 1)
            age = sim.t-e0
            born = pas2[hood].any(1)
            kept = ~born & opening[hood].any(1) & (age < RIFT_STALL_MYR)
            hk = hood[kept]
            covered[hk[opening[hk]]] = True
            b['end_passive_km'] += float(born.sum())*self.sp_km
            b['passive_age_sum_km_myr'] += float(age[born].sum())*self.sp_km
            dead = ~born & ~kept
            jd = j[dead]; hd = hood[dead]
            is_edge = edge[hd].any(1); is_cont = cont[hd].any(1)
            fate = np.full(jd.size, 'other', dtype=object)
            fate[is_cont & is_edge & near[jd] & (nearest[jd] == SUBDUCTION)] = 'subduction'
            fate[is_cont & is_edge & near[jd] & (nearest[jd] == COLLISION)] = 'collision'
            fate[is_cont & ~is_edge] = 'closed'
            fate[~is_cont] = 'lost'
            fate[opening[hd].any(1) & (age[dead] >= RIFT_STALL_MYR)] = 'stalled'
            for f in RIFT_FATES[1:]:
                b[f'end_{f}_km'] += float((fate == f).sum())*self.sp_km
            xs.append(x[kept]); es.append(e0[kept]); js.append(j[kept])
        oage = np.where(cont, np.inf, np.nan_to_num(np.asarray(st['ocean_age_myr'], np.float64), nan=np.inf))
        young1 = np.minimum(oage, np.where(ok, oage[nbc], np.inf).min(1))
        young = np.minimum(young1, np.where(ok, young1[nbc], np.inf).min(1)) < RIFT_YOUNG_MYR
        new = np.flatnonzero(edge & near & (nearest == RIDGE) & young & ~covered)
        b['entered_km'] += float(new.size)*self.sp_km
        xs.append(g.xyz[new]); es.append(np.full(new.size, float(sim.t))); js.append(new)
        x = np.concatenate(xs); e = np.concatenate(es); jn = np.concatenate(js)
        b['tracked_km'] += float(x.shape[0])*self.sp_km
        rot = {pl.id: fields.rotation(pl.axis, pl.speed*1e-2*sim.p['dt_myr']*1e6/g.radius) for pl in sim.plates}
        own = sim.owner[jn]
        for pid in np.unique(own):
            if pid in rot:
                m = own == pid
                x[m] = x[m]@rot[pid].T
        self.prev = (x, e)

    def summary(self):
        return {'bin_myr': BIN_MYR, 'node_km': self.sp_km, 'young_myr': RIFT_YOUNG_MYR, 'stall_myr': RIFT_STALL_MYR,
                'bins': [self.bins[k] for k in sorted(self.bins)]}


def pool_rift(records, t_from=-1e9, t_to=1e9):
    tot = {}; steps = 0
    for _, rec in records:
        for b in rec['bins']:
            if t_from <= b['t0_myr'] < t_to:
                steps += b['steps']
                for k, v in b.items():
                    if k not in ('t0_myr', 'steps'):
                        tot[k] = tot.get(k, 0.)+v
    if not steps:
        return None
    n = len(records); myr = steps/n
    ended = sum(tot[f'end_{f}_km'] for f in RIFT_FATES)
    return {'myr': myr, 'entered_km_myr': tot['entered_km']/n/myr, 'tracked_km': tot['tracked_km']/steps,
            'ended_km_myr': ended/n/myr, 'conversion': tot['end_passive_km']/max(ended, 1e-9),
            'fates': {f: tot[f'end_{f}_km']/max(ended, 1e-9) for f in RIFT_FATES},
            'time_to_passive_myr': tot['passive_age_sum_km_myr']/max(tot['end_passive_km'], 1e-9)}


def fmt_rift(name, r):
    if r is None:
        return f'{name}: no rift record'
    return (f"{name}: young rifted margin enters {r['entered_km_myr']:5.0f} km/Myr (tracked {r['tracked_km']/1e3:5.1f}k km), "
            f"ends {r['ended_km_myr']:5.0f} km/Myr; becomes passive {r['conversion']:.2f} of it (Earth ~0.6-0.85, steady "
            f"state), after {r['time_to_passive_myr']:4.0f} Myr\n      fates: " + ', '.join(f'{k} {v:.2f}' for k, v in r['fates'].items()))


def run_dirs(paths):
    """Seed folders (with stats.json) under the given sweep, combination or seed folders."""
    out = []
    for f in map(Path, paths):
        out += [f] if (f/'stats.json').is_file() else sorted(p for p in f.rglob('seed*') if (p/'stats.json').is_file())
    return out


def read_runs(paths, key='passive_record'):
    groups = {}
    for d in run_dirs(paths):
        st = json.loads((d/'stats.json').read_text())
        if key in st:
            groups.setdefault(d.parent, []).append((d.name, st[key]))
    return groups


def pool(records, t_from=-1e9, t_to=1e9):
    """Sums over the bins in [t_from, t_to) and the seeds; passive km as the mean over steps."""
    tot = {}
    steps = 0
    for _, rec in records:
        for b in rec['bins']:
            if not (t_from <= b['t0_myr'] < t_to):
                continue
            steps += b['steps']
            for k, v in b.items():
                if k in ('t0_myr', 'steps'):
                    continue
                tot[k] = (np.add(tot[k], v) if isinstance(v, list) else tot.get(k, 0.)+v) if k in tot else (list(v) if isinstance(v, list) else v)
    n = len(records)
    if not steps:
        return None
    myr = steps/n
    ended = sum(tot[f'end_{f}_km'] for f in FATES)
    out = {'myr': myr, 'passive_km': tot['passive_km']/steps,
           'births_km_myr': (tot['births_rifted_km']+tot['births_converted_km'])/n/myr,
           'births_rifted_share': tot['births_rifted_km']/max(tot['births_rifted_km']+tot['births_converted_km'], 1e-9),
           'ended_km_myr': ended/n/myr,
           'fates': {f: tot[f'end_{f}_km']/max(ended, 1e-9) for f in FATES},
           # records written before the hysteresis (6 Nov pre-check) have no lasting ends: zeros
           'lasting_fates': {f: tot.get(f'end_{f}_lasting_km', 0.)/max(sum(tot.get(f'end_{g}_lasting_km', 0.) for g in FATES), 1e-9) for f in FATES},
           'end_by_age_share': [v/max(sum(tot['end_by_age_km']), 1e-9) for v in tot['end_by_age_km']],
           'durable_end_km_myr': sum(tot['end_by_age_km'][1:])/n/myr,
           'mean_age_at_end_myr': tot['end_age_sum_km_myr']/max(sum(tot['end_by_age_km']), 1e-9),
           'start_cohort_end_km_myr': tot['end_start_cohort_km']/n/myr}
    out['tracked_km'] = tot.get('tracked_km', tot['passive_km'])/steps
    out['lifespan_myr'] = out['tracked_km']/max(out['ended_km_myr'], 1e-9)
    out['durable_lifespan_myr'] = out['tracked_km']/max(out['durable_end_km_myr'], 1e-9)
    return out


def fmt(name, r):
    if r is None:
        return f'{name}: no passive record'
    f = r['fates']
    return (f"{name}: passive {r['passive_km']/1e3:5.1f}k km, births {r['births_km_myr']:5.0f} km/Myr (rifted {r['births_rifted_share']:.2f}), "
            f"ends {r['ended_km_myr']:5.0f} km/Myr (lasting >= 5 Myr {r['durable_end_km_myr']:5.0f}), lifespan "
            f"{r['lifespan_myr']:4.0f} Myr (lasting ends {r['durable_lifespan_myr']:4.0f}; Earth 130-181), mean age at end "
            f"{r['mean_age_at_end_myr']:4.0f} Myr\n      ends: " + ', '.join(f'{k} {v:.2f}' for k, v in f.items())
            + '\n      ends lasting >= 5 Myr (and the start cohort): ' + ', '.join(f'{k} {v:.2f}' for k, v in r['lasting_fates'].items())
            + '\n      age at end ' + ' / '.join(f'{a:.0f}-{b:.0f}' if np.isfinite(b) else f'{a:.0f}+' for a, b in zip(AGE_EDGES[:-1], AGE_EDGES[1:]))
            + ' Myr: ' + ' / '.join(f'{v:.2f}' for v in r['end_by_age_share']))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--json')
    a = ap.parse_args(argv)
    print(f"Earth (Bradley 2008): ancient passive margins' mean lifespan 181 Myr (Permian-Neogene 130); today's mean age "
          f"104 Myr; passive edge {EARTH_PASSIVE_KM/1e3:.1f}k km at 46 km -> ~{EARTH_PASSIVE_KM/181:.0f}-{EARTH_PASSIVE_KM/130:.0f} km/Myr ending")
    out = {}
    for grp, recs in read_runs(a.paths).items():
        name = f'{grp.name} ({len(recs)} seeds)'
        whole, late = pool(recs), pool(recs, -150., 1e9)
        print(fmt(name+', whole run', whole)); print(fmt('   last 150 Myr', late))
        for s, rec in recs:
            r = pool([(s, rec)], -150., 1e9)
            print(f"   {s}: last 150 Myr passive {r['passive_km']/1e3:.1f}k km, ends {r['ended_km_myr']:.0f} km/Myr, lifespan "
                  f"{r['lifespan_myr']:.0f} Myr, subduction {r['fates']['subduction']:.2f}, collision+closed "
                  f"{r['fates']['collision']+r['fates']['closed_collision']+r['fates']['closed']:.2f}, ridge {r['fates']['ridge']:.2f}, lost {r['fates']['lost']:.2f}")
        out[str(grp)] = {'whole': whole, 'last150': late}
    for grp, recs in read_runs(a.paths, 'rift_record').items():
        name = f'{grp.name} ({len(recs)} seeds)'
        whole, late = pool_rift(recs), pool_rift(recs, -150., 1e9)
        print(fmt_rift('rift record '+name+', whole run', whole)); print(fmt_rift('   last 150 Myr', late))
        for s, rec in recs:
            r = pool_rift([(s, rec)], -150., 1e9)
            print(f"   {s}: last 150 Myr enters {r['entered_km_myr']:.0f} km/Myr, passive {r['conversion']:.2f}, subduction "
                  f"{r['fates']['subduction']:.2f}, collision+closed {r['fates']['collision']+r['fates']['closed']:.2f}, stalled {r['fates']['stalled']:.2f}")
        out.setdefault(str(grp), {})['rift'] = {'whole': whole, 'last150': late}
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
