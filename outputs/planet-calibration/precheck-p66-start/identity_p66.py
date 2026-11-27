import sys, json, hashlib, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
params = json.loads(sys.argv[2]); record = len(sys.argv) > 3 and sys.argv[3] == 'record'
g = CubeSphere(48)
out = []
for seed in (3, 11):
    sim = PlateSim(g, seed=seed, params={'threads': 1, **params})
    st0 = sim.state()                              # as simulate_plates' first record (state() sets z_last)
    if record:
        from planet_reviews.passive_survival import PassiveRecorder
        rec = PassiveRecorder(g); rec.update(sim, st0)
    for s in range(70):
        sim.step()
        st = sim.state()
        if record:
            rec.update(sim, st)
        h = hashlib.md5()
        for k in sorted(st):
            v = np.asarray(st[k])
            if v.dtype.kind in 'fiub':
                h.update(np.ascontiguousarray(v).tobytes())
        out.append(h.hexdigest())
open(sys.argv[1], 'w').write('\n'.join(out))
print(len(out), 'steps', 'record' if record else '')
