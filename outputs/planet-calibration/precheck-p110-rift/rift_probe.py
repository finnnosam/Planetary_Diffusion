"""Pass 110 (journal 148): the rift's volume balance, read per step from a warm state (STATE section 0 item 1; V3).

Question (fixed before the runs): does the crust the new plate nodes receive at continental rift gaps fall short of the
volume the world's pure shear takes from the stretching zone, by about the shift of V3's update part between 'nearest'
and 'overlap' (-0.27 km^3/yr on the Remap Sweep's six seeds)?

Per step (hooks on `_rift_pure_shear` and `_update_rasters`; the trajectory is not changed):
  thinned      the zone's pure-shear thinning on the world, -sum a dT (m^3)
  rift_world   the gap crust on the world, sum_r Tg_r a_r (equal to thinned: the kernel conserves on the world grid)
  A/a          the area the plates open over each rift gap (unmasked plate nodes of the owning plate whose forward
               nearest node is the gap) against the gap's own area: mean, share of gaps with A = 0, with A > a
  placed       what the new plate nodes over rift gaps receive: sum_r A_r T_end,r (T_end the gap's end-of-step thickness)
  placed_rift  the rift part of it, sum_r A_r Tg_r; shortfall = placed_rift - rift_world
  other_gap    the rest at the gap, sum_r (A_r - a_r)(T_end,r - Tg_r) (flow, erosion at gap nodes: unbooked by nearest)
  new_ledger   the series' new-node ledger step (continental), for the check placed ~ new_ledger

  python rift_probe.py <warm dir> <remap_rule> [steps] [--param k=v ...]
"""
import copyreg
import importlib
import json
import pickle
import sys
import types

import numpy as np

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))
from planet_grid import fields  # noqa: E402
from planet_tectonics.params import PARAMS  # noqa: E402

args = [x for x in sys.argv[1:] if not x.startswith('--param')]
pars = [x.split('=', 1) for i, x in enumerate(sys.argv) if i > 0 and sys.argv[i-1] == '--param']
args = [x for x in args if '=' not in x]
warm, rule = args[0], args[1]
K = int(args[2]) if len(args) > 2 else 6
sim = pickle.load(open(f'{warm}/warm.pkl', 'rb'))
for k, v in PARAMS.items():
    sim.p.setdefault(k, v)
sim.p['remap_rule'] = rule
for k, v in pars:
    sim.p[k] = json.loads(v) if v[:1] in '0123456789-[{tf"' and v not in ('plates', 'world') else v
sim.state()
g = sim.g
a = g.area
rec = {}
o_rps, o_up = sim._rift_pure_shear, sim._update_rasters


def rps(rows, ncont, solid, cn, nbs):
    w = sim.world
    c0 = w['cont'].copy()
    th0 = w['thick'].astype(np.float64).copy()
    dT = o_rps(rows, ncont, solid, cn, nbs)
    r = np.flatnonzero(w['cont'] & ~c0 & sim.new_crust)
    rec['rows'] = r
    rec['Tg'] = w['thick'][r].astype(np.float64).copy()
    rec['thinned'] = -float((a*dT).sum())
    rec['rift_world'] = float((a[r]*rec['Tg']).sum())
    return dT


def up(dthick, idx, dsed=None, drock=None):
    w = sim.world
    A = np.zeros(g.size)
    for k, pl in enumerate(sim.plates):
        wn = fields.nearest_nodes(g, g.xyz@pl.Q.T)
        new = (sim.owner_slot[wn] == k) & ~pl.mask & sim.new_crust[wn]
        np.add.at(A, wn[new], a[new])
    r = rec.get('rows', np.zeros(0, np.int64))
    Tg = rec.get('Tg', np.zeros(0))
    Te = w['thick'][r].astype(np.float64)
    Ar, ar = A[r], a[r]
    rec['n_rift'] = int(r.size)
    rec['A_over_a'] = float(Ar.sum()/max(ar.sum(), 1e-9))
    rec['zeroA'] = float((Ar == 0).mean()) if r.size else 0.
    rec['bigA'] = float((Ar > 1.01*ar).mean()) if r.size else 0.
    rec['placed'] = float((Ar*Te).sum())
    rec['placed_rift'] = float((Ar*Tg).sum())
    rec['shortfall'] = rec['placed_rift']-rec.get('rift_world', 0.)
    rec['other_gap'] = float(((Ar-ar)*(Te-Tg)).sum())
    L0 = sim.ledger.get('new_nodes', 0.)
    out = o_up(dthick, idx, dsed, drock)
    rec['new_ledger'] = sim.ledger.get('new_nodes', 0.)-L0
    return out


sim._rift_pure_shear, sim._update_rasters = rps, up
tot = {}
for i in range(K):
    rec.clear()
    sim.step()
    sim.state()
    row = {k: v for k, v in rec.items() if k not in ('rows', 'Tg')}
    for k in ('thinned', 'rift_world', 'placed', 'placed_rift', 'shortfall', 'other_gap', 'new_ledger'):
        row[k] = row.get(k, 0.)/1e15                      # 1e6 km^3
        tot[k] = tot.get(k, 0.)+row[k]
    print(f"{rule} t={sim.t:+.0f}: " + ', '.join(f'{k} {v:+.4g}' if isinstance(v, float) else f'{k} {v}'
                                             for k, v in row.items()), flush=True)
print(f"{rule} total over {K} steps (1e6 km^3): " + ', '.join(f'{k} {v:+.4f}' for k, v in tot.items()))
print(f"{rule} shortfall per Myr {tot['shortfall']/K:+.3f} km^3/yr; other_gap {tot['other_gap']/K:+.3f}; "
      f"placed - new_ledger {(tot['placed']-tot['new_ledger'])/K:+.3f}")
