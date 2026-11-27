import sys, json, hashlib, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
from planet_reviews.passive_survival import PassiveRecorder
params = json.loads(sys.argv[2])
g = CubeSphere(48)
out = []
for seed in (3, 11):
    sim = PlateSim(g, seed=seed, params={'threads': 1, **params})
    st0 = sim.state(); rec = PassiveRecorder(g); rec.update(sim, st0)
    for s in range(70):
        sim.step(); st = sim.state(); rec.update(sim, st)
        h = hashlib.md5()
        for k in sorted(st):
            v = np.asarray(st[k])
            if v.dtype.kind in 'fiub':
                h.update(np.ascontiguousarray(v).tobytes())
        out.append(h.hexdigest())
    out.append('events ' + json.dumps([(e['t'], e['kind']) for e in sim.events][-50:] if hasattr(sim, 'events') else []))
open(sys.argv[1], 'w').write('\n'.join(out))
print(len(out), 'lines')
