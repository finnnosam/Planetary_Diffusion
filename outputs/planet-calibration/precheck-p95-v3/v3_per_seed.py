"""V3 per seed: crust_budget's closure residuals (crust; ocean-floor sediment), whole run and last 40 Myr (pass 95)."""
import json, sys
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import crust_budget as cb
g = CubeSphere(217)
for arm in sys.argv[1:]:
    print(Path(arm).name)
    for s in sorted(Path(arm).glob('seed*')):
        st = json.loads((s/'stats.json').read_text()); p = cb.run_params(s)
        T0 = st['series'][0]['time_myr']
        w, r = cb.budget(s, g, T0, 0., p, st), cb.budget(s, g, -40., 0., p, st)
        print(f"  {s.name}: whole crust {w['crust']['residual']:+.3f} floor {w['ocean_floor_sediment']['residual']:+.3f} | "
              f"last 40 crust {r['crust']['residual']:+.3f} floor {r['ocean_floor_sediment']['residual']:+.3f} km3/yr")
