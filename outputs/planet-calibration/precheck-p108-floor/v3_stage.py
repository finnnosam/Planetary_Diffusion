"""Pass 108 (journal 145): one step from the warm state, the world's continental volume by stage (composition, after
rifting / flow / thickening, after erosion), dthick's volume, and the plates' change by the increments, per rule.

  python v3_stage.py <warm dir> <rule> [steps]
"""
import copyreg
import importlib
import pickle
import sys
import types

import numpy as np

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))
from planet_tectonics.params import PARAMS  # noqa: E402

warm, rule = sys.argv[1], sys.argv[2]
K = int(sys.argv[3]) if len(sys.argv) > 3 else 2
sim = pickle.load(open(f'{warm}/warm.pkl', 'rb'))
for k, v in PARAMS.items():
    sim.p.setdefault(k, v)
sim.p['remap_rule'] = rule
sim.state()
g = sim.g
a = g.area
rec = {}
o_wf, o_er, o_up, o_rift = sim._world_from, sim._erode, sim._update_rasters, sim._rift


def cvol(th, cont):
    return float((a*th*cont).sum())/1e9


def plate_cont():
    return sum(float((a*pl.thick.astype(np.float64)*(pl.mask & pl.cont)).sum()) for pl in sim.plates)/1e9


def wf(idx, claimed, record=True):
    out = o_wf(idx, claimed, record)
    if record:
        w = sim.world
        rec['cont0'] = w['cont'].copy()
        rec['V_composed'] = cvol(w['thick'], w['cont'] & ~sim.new_crust)
        rec['V_plate0'] = plate_cont()
    return out


def rift():
    w = sim.world
    v0 = cvol(w['thick'], w['cont'])
    out = o_rift()
    rec['rift_dthick_vol'] = float((a*np.asarray(out)*rec['cont0']).sum())/1e9 if out is not None else 0.
    return out


def er(th):
    w = sim.world
    rec['V_tect'] = cvol(th, rec['cont0'] & ~sim.new_crust)
    e, ds = o_er(th)
    rec['erosion_vol'] = float((a*e*rec['cont0']).sum())/1e9
    return e, ds


def up(dthick, idx, dsed, drock):
    w = sim.world
    c = w['cont'] & ~sim.new_crust
    rec['dthick_cont_now'] = float((a*dthick*c).sum())/1e9
    rec['dthick_cont0'] = float((a*dthick*rec['cont0'] & ~sim.new_crust).sum())/1e9 if False else float((a*dthick*(rec['cont0'] & ~sim.new_crust)).sum())/1e9
    rec['cont_flips'] = int((rec['cont0'] != w['cont']).sum())
    L0 = sim.ledger.get('increments', 0.)
    out = o_up(dthick, idx, dsed, drock)
    rec['ledger_increments'] = (sim.ledger.get('increments', 0.)-L0)/1e9
    return out


sim._world_from, sim._erode, sim._update_rasters, sim._rift = wf, er, up, rift
for i in range(K):
    E0 = sim.erosion_totals.get('continental_crust_eroded_m3', 0.)
    rec.clear()
    sim.step()
    sim.state()
    rec['booked_erosion'] = (sim.erosion_totals.get('continental_crust_eroded_m3', 0.)-E0)/1e9
    print(f"{rule} t={sim.t:+.0f}: " + ', '.join(f'{k} {v:+.4g}' if isinstance(v, float) else f'{k} {v}'
                                             for k, v in rec.items() if k not in ('cont0',)), flush=True)
