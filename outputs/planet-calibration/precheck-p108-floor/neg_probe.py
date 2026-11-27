import copyreg, importlib, pickle, sys, types
import numpy as np
copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))
from planet_tectonics.params import PARAMS
from planet_grid import fields
sim = pickle.load(open('/home/claude/t1028/warm.pkl','rb'))
for k, v in PARAMS.items(): sim.p.setdefault(k, v)
sim.p['remap_rule']='overlap'
sim.state()
g=sim.g
track=None
for i in range(8):
    sim.step(); sim.state()
    rows=[]
    for k,pl in enumerate(sim.plates):
        wn=fields.nearest_nodes(g, g.xyz@pl.Q.T)
        own=sim.owner_slot[wn]==k
        neg=np.flatnonzero(pl.mask & (pl.thick<0))
        for j in neg: rows.append((float(pl.thick[j]), k, int(j), bool(pl.cont[j]), bool(own[j]), float(sim.world['thick'][wn[j]]), int(sim.owner_slot[wn[j]]), float(sim.t-pl.born[j]) if np.isfinite(pl.born[j]) else -1))
    rows.sort()
    print('t=%+.0f negative plate nodes %d; cont %d; owned %d; worst:'%(sim.t, len(rows), sum(r[3] for r in rows), sum(r[4] for r in rows)), [(round(r[0]),r[3],r[4],round(r[5]),round(r[7],1)) for r in rows[:4]], flush=True)
    print('   world thick min %.0f; world cont thick min %.0f'%(sim.world['thick'].min(), sim.world['thick'][sim.world['cont']].min()), flush=True)
