"""Pass 108 (journal 145): where the remap's crust-budget residual comes from. From the warm state, K steps per rule;
per step the plates' continental crust ledger by key and crust_budget's booked terms (erosion_totals), km^3 per step.

  python v3_probe.py <warm dir> <rule> [K]
"""
import copyreg
import importlib
import json
import pickle
import sys
import types

import numpy as np

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))
from planet_tectonics.params import PARAMS  # noqa: E402

TERMS = ('continental_crust_eroded_m3', 'subduction_on_continents_m3', 'sediment_accreted_m3',
         'collision_sediment_to_crust_m3', 'volcanic_added_m3', 'arc_accreted_on_continent_m3', 'eclogite_foundered_m3')

warm, rule = sys.argv[1], sys.argv[2]
K = int(sys.argv[3]) if len(sys.argv) > 3 else 4
sim = pickle.load(open(f'{warm}/warm.pkl', 'rb'))
for k, v in PARAMS.items():
    sim.p.setdefault(k, v)
sim.p['remap_rule'] = rule
sim.state()
g = sim.g
A = 4*np.pi*(g.radius/1e3)**2
rows = []
for i in range(K):
    L0, E0, V0 = dict(sim.ledger), dict(sim.erosion_totals), sim.continental_volume()
    AL0 = dict(sim.area_ledger)
    sim.step()
    sim.state()
    dL = {k: (sim.ledger.get(k, 0.)-L0.get(k, 0.))/1e9 for k in set(sim.ledger) | set(L0)}
    dE = {k: (sim.erosion_totals.get(k, 0.)-E0.get(k, 0.))/1e9 for k in TERMS}
    arc = (sim.area_ledger.get('arc', 0.)-AL0.get('arc', 0.))/1e6*sim.p['arc_to_continent_thick_m']/1e3   # m^2 -> km^3 (crust_budget's lower bound)
    booked = sum(dE.values())-2*dE['eclogite_foundered_m3']+arc
    net = (sim.continental_volume()-V0)/1e9
    rows.append({'t': sim.t, 'net': net, 'booked': booked, 'residual': net-booked, 'ledger': dL, 'terms': dE, 'arc_conv': arc})
    print(f"{rule} t={sim.t:+.0f}: net {net:+.1f} km3, booked {booked:+.1f}, residual {net-booked:+.1f}; ledger "
          + ', '.join(f'{k} {v:+.1f}' for k, v in sorted(dL.items()) if abs(v) > .05)
          + '; terms ' + ', '.join(f'{k.replace("_m3", "")} {v:+.1f}' for k, v in dE.items() if abs(v) > .05) + f'; arc_conv {arc:+.1f}', flush=True)
json.dump(rows, open(f'{warm}/v3_{rule}.json', 'w'), indent=1)
