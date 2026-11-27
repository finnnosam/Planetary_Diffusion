"""Pass 103 probe: a full run with every split_rule 'earth' draw logged (no change to the model: the same code, the same
rng draws, in the same order). Per plate with continent and per step: t, plate, plate area, continental area Ac, the
line's kind (force / chord / bisector), its margin L (km), prod, the acceptance probability, accepted, split made, the
new part's area, the candidates tried, and the plate's boundary make-up (trench / overriding / ridge nodes on its edge).
Usage (repo root): python3 split_log_run.py <out> <seed> [key=value ...]"""
import json, sys
import numpy as np
from planet_tectonics import events as ev
from planet_tectonics.params import RIDGE, SUBDUCTION

LOG = []


def _split_events_earth(self):
    from planet_tectonics.breakup import production_per_mkm2
    g, p, rng, w = self.g, self.p, self.rng, self.world
    dt = p['dt_myr']
    A = g.area
    prod = production_per_mkm2(self.t-self.t_start)
    for pl in list(self.plates):
        own = self.owner == pl.id
        area = A[own].sum()/1e12
        if area < p['split_min_area_mkm2']:
            continue
        Ac = A[own & w['cont']].sum()/1e12
        done = False
        rec = None
        if prod*Ac > 0:
            lab = self._rift_line(pl, own)
            L = self._split_margin_km(own, lab)
            pr = 1-np.exp(-prod*Ac/max(L, ev.SPLIT_L_MIN_KM)*dt) if L > 0 else 0.
            edge = own & np.any((g.nbr >= 0) & (self.owner[np.clip(g.nbr, 0, None)] != pl.id), 1)
            rec = {'t': float(self.t), 'plate': int(pl.id), 'area': float(area), 'Ac': float(Ac), 'prod': float(prod),
                   'kind': self._line_kind, 'tries': int(getattr(self, '_chord_tries', -1)), 'L': float(L), 'p': float(pr),
                   'part': float(A[own & (lab == 1)].sum()/1e12), 'part_c': float(A[own & (lab == 1) & w['cont']].sum()/1e12),
                   'edge': int(edge.sum()), 'trench': int((edge & self.trench).sum()),
                   'over': int((edge & (self.btype == SUBDUCTION) & ~self.trench).sum()),
                   'ridge': int((edge & (self.btype == RIDGE)).sum()), 'acc': False, 'made': False}
            if L > 0 and rng.uniform() < 1-np.exp(-prod*Ac/max(L, ev.SPLIT_L_MIN_KM)*dt):
                rec['acc'] = True
                if self.split(pl, lab == 1, self._rift_speed(L), 'rift') is not None:
                    rec['made'] = True
                    self.rift_margin_km = getattr(self, 'rift_margin_km', 0.)+L
                    if p['rift_line'] in ('chord', 'force'):
                        self.events[-1]['line'] = self._line_kind
                done = True
            LOG.append(rec)
        if not done:
            rate = p['split_rate_per_myr']*area/100.
            if rng.uniform() < 1-np.exp(-rate*dt):
                lab = self._rift_line(pl, own, chord=False)
                if self._split_margin_km(own, lab) == 0:
                    self.split(pl, lab == 1, p['split_speed_cm_yr'], 'rift')
                    done = True
        if done:
            continue
        old = own & ~w['cont'] & (self.t-w['born'] > p['init_age_myr'])
        old_area = A[old].sum()/1e12
        irate = p['init_rate_per_myr']*min(5., old_area/p['init_area_ref_mkm2'])
        if old.sum() >= 4 and rng.uniform() < 1-np.exp(-irate*dt):
            self._initiate_subduction(pl, own, old)


ev.EventsMixin._split_events_earth = _split_events_earth
if __name__ == '__main__':
    out, seed = sys.argv[1], sys.argv[2]
    sys.path.insert(0, 'scripts')
    import simulate_plates as sp
    args = ['--seeds', seed, '--out', out, '--save-every', '5']
    for kv in sys.argv[3:]:
        args += ['--param', kv]
    try:
        sp.main(args)
    finally:
        from pathlib import Path
        Path(out).mkdir(parents=True, exist_ok=True)
        (Path(out)/f'splitlog{seed}.json').write_text(json.dumps(LOG))
