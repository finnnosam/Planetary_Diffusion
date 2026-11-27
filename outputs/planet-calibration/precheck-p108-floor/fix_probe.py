import copyreg, importlib, pickle, sys, types
import numpy as np
copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))
from planet_tectonics.params import PARAMS
rule=sys.argv[1]
sim = pickle.load(open('/home/claude/t1028/warm.pkl','rb'))
for k, v in PARAMS.items(): sim.p.setdefault(k, v)
sim.p['remap_rule']=rule
sim.state()
def rng():
    th=np.concatenate([pl.thick[pl.mask].astype(float) for pl in sim.plates]); sd=np.concatenate([pl.sed[pl.mask].astype(float) for pl in sim.plates])
    return th.min(), np.percentile(th,99.99), th.max(), sd.max()
print(rule, 'start: thick min %.0f p99.99 %.0f max %.0f, sed max %.0f'%rng())
for i in range(8):
    E0=dict(sim.erosion_totals); sim.step(); sim.state()
    mv={k:(sim.erosion_totals.get(k,0)-E0.get(k,0))/1e9 for k in sim.erosion_totals if k.startswith('return_rescaled')}
    print(rule, 'clamped km3 %.1f'%(sim.erosion_totals.get('return_thick_clamped_m3',0)/1e9), 't=%+.0f thick min %.0f p99.99 %.0f max %.0f, sed max %.0f'%((sim.t,)+rng()), 'moved km3', {k[16:]:round(v,1) for k,v in mv.items()}, 'cont vol %.6g'%sim.continental_volume(), flush=True)
