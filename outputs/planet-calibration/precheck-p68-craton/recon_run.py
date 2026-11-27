"""Pass 68 probe (cloud only, output only): reconcile the rift record ('becomes passive') with the survival record's
births. Runs scripts.simulate_plates with a patched copy of planet_reviews.passive_survival that also records, each
Myr: the rift particles that become passive (km), the distinct passive-edge nodes within their reach (two rings) (km),
how many of those are the survival record's births this Myr, and how many were passive or kept already; and the
survival births not reached by any newly passive particle. Writes <out>/seed<seed>/recon.json.
usage: python recon_run.py <seed> <out> [--param k=v ...]"""
import json, sys, types
import numpy as np
import planet_reviews
from pathlib import Path

src = Path(planet_reviews.__file__).with_name('passive_survival.py').read_text()
a1 = "        birth[new] = sim.t\n"
assert src.count(a1) == 1
src = src.replace(a1, a1+"        self.new_last = new.copy(); self.reached_last = reached.copy()\n")
a2 = "            b['end_passive_km'] += float(born.sum())*self.sp_km\n"
assert src.count(a2) == 1
src = src.replace(a2, a2+"""            D = getattr(self, 'diag', None)
            if D is None:
                D = self.diag = {}
            hb = hood[born]
            reach = np.zeros(g.size, bool); reach[hb.ravel()] = True
            reach2 = reach | np.any(np.where(ok, reach[nbc], False), 1)
            hit = reach2 & pas
            newp = getattr(self, '_passive_new', np.zeros(g.size, bool))
            d = D.setdefault(float(np.floor(sim.t/10.)*10.), {'born_km': 0., 'hit_km': 0., 'hit_new_km': 0., 'hit_old_km': 0.,
                                                              'births_km': 0., 'births_unreached_km': 0., 'passive_km': 0., 'steps': 0})
            d['steps'] += 1
            d['born_km'] += float(born.sum())*self.sp_km
            d['hit_km'] += float(hit.sum())*self.sp_km
            d['hit_new_km'] += float((hit & newp).sum())*self.sp_km
            d['hit_old_km'] += float((hit & ~newp).sum())*self.sp_km
            d['births_km'] += float(newp.sum())*self.sp_km
            d['births_unreached_km'] += float((newp & ~reach2).sum())*self.sp_km
            d['passive_km'] += float(pas.sum())*self.sp_km
""")
mod = types.ModuleType('planet_reviews.passive_survival')
mod.__file__ = planet_reviews.__file__
exec(compile(src, 'passive_survival_probe', 'exec'), mod.__dict__)
sys.modules['planet_reviews.passive_survival'] = mod
planet_reviews.passive_survival = mod

# hand the passive recorder's births of this Myr to the rift recorder (it runs right after it, same state)
_upd = mod.RiftRecorder.update
holder = {}
def upd(self, sim, st, reading):
    self._passive_new = holder['p'].new_last if 'p' in holder and hasattr(holder['p'], 'new_last') else np.zeros(self.g.size, bool)
    holder['r'] = self
    return _upd(self, sim, st, reading)
mod.RiftRecorder.update = upd
_pinit = mod.PassiveRecorder.__init__
def pinit(self, g):
    _pinit(self, g); holder['p'] = self
mod.PassiveRecorder.__init__ = pinit

seed, out = sys.argv[1], sys.argv[2]
from scripts import simulate_plates
simulate_plates.main(['--seeds', seed, '--out', out, '--save-every', '5', '--param', 'threads=1']+sys.argv[3:])
Path(out, f'seed{seed}', 'recon.json').write_text(json.dumps(getattr(holder['r'], 'diag', {}), indent=1))
