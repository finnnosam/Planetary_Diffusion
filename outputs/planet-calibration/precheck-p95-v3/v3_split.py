"""V3 split (pass 95): crust_budget's continental-crust residual per seed, split by the plate layer's own crust ledger,
which runs since pass 95 write to the stats series ('crust_ledger_km3'; sim.ledger, whose sum is the crust volume's
change). For the whole run and the last 40 Myr, km^3/yr:

  closure    = d(volume) - d(sum of the ledger)          (the ledger's own closure: ~0)
  residual   = d(volume) - crust_budget's booked terms  (as crust_budget prints it)
             = merge + pieces + holes + update + closure
  update     = the rest: what the raster update applied (increments, ownership handovers: new nodes, collision
               transfers, back-arc strips, arcs) against the booked terms (erosion, trench, sediment accreted,
               collision sediment, volcanic, arcs, foundering); its gross parts run at several km^3/yr and cancel, so
               only its net is read, also as a share of crystalline erosion (the nearest-node writeback's loss).

Runs without the key print 'no crust ledger'. Usage: v3_split.py <arm folders>."""
import json
import sys
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_reviews import crust_budget as cb

EVENTS = ('merge', 'pieces', 'holes')
g = CubeSphere(217)


def split(seed, t0, t1):
    st = json.loads((seed/'stats.json').read_text()); p = cb.run_params(seed)
    S = {round(r['time_myr'], 3): r for r in st['series']}
    a, b = S[round(t0, 3)], S[round(t1, 3)]
    if 'crust_ledger_km3' not in b:
        return None
    yr = (t1-t0)*1e6
    la, lb = a.get('crust_ledger_km3', {}), b['crust_ledger_km3']
    dl = {k: (lb.get(k, 0.)-la.get(k, 0.))/yr for k in set(la) | set(lb)}
    dv = (b['continental_volume_km3']-a['continental_volume_km3'])/yr
    bud = cb.budget(seed, g, t0, t1, p, st)['crust']
    res = bud['residual']
    ev = {k: dl.get(k, 0.) for k in EVENTS}
    closure = dv-sum(dl.values())
    wb = res-sum(ev.values())-closure
    return {'residual': res, **ev, 'update': wb, 'closure': closure, 'erosion': bud['crystalline_erosion'],
            'wb_share': wb/abs(bud['crystalline_erosion']) if bud['crystalline_erosion'] else np.nan}


for arm in sys.argv[1:]:
    print(Path(arm).name)
    for seed in sorted(Path(arm).glob('seed*')):
        st = json.loads((seed/'stats.json').read_text()); T0 = st['series'][0]['time_myr']
        for nm, (t0, t1) in (('whole run', (T0, 0.)), ('last 40', (max(-40., T0), 0.))):
            r = split(seed, t0, t1)
            if r is None:
                print(f'  {seed.name} {nm}: no crust ledger in the series (run before pass 95)')
                continue
            print(f"  {seed.name} {nm:9s}: residual {r['residual']:+.3f} = merge {r['merge']:+.3f} + pieces {r['pieces']:+.3f} + "
                  f"holes {r['holes']:+.3f} + update {r['update']:+.3f} + closure "
                  f"{r['closure']:+.4f} km3/yr; update / crystalline erosion {r['wb_share']:+.3f}")
