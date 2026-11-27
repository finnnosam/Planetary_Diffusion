"""Per-step crust residual split by the step's events (merge / detach / rift / initiation counts rising), pass 95."""
import json, sys
from pathlib import Path
import numpy as np
from planet_reviews.crust_review import run_params
R = 6371.0
KEYS = ['continental_crust_eroded_m3', 'subduction_on_continents_m3', 'sediment_accreted_m3',
        'collision_sediment_to_crust_m3', 'volcanic_added_m3', 'arc_accreted_on_continent_m3']
for arm in sys.argv[1:]:
    rows = []
    for seed in sorted(Path(arm).glob('seed*')):
        st = json.loads((seed/'stats.json').read_text()); p = run_params(seed); A = 4*np.pi*R**2
        S = st['series']
        for a, b in zip(S[:-1], S[1:]):
            E = lambda k: b['erosion_totals_km3'].get(k, 0.)-a['erosion_totals_km3'].get(k, 0.)
            arc = (b['area_ledger_frac'].get('arc', 0.)-a['area_ledger_frac'].get('arc', 0.))*A*p['arc_to_continent_thick_m']/1e3
            r = b['continental_volume_km3']-a['continental_volume_km3']-(sum(E(k) for k in KEYS)+arc-E('eclogite_foundered_m3'))
            ev = {k: b['events'].get(k, 0)-a['events'].get(k, 0) for k in ('merge', 'detach', 'rift', 'subduction initiation')}
            rows.append((r/1e6, ev))
    r = np.array([x[0] for x in rows])
    big = np.abs(r) > 3
    print(f"{Path(arm).name}: steps {len(r)}, |res| > 3 Mkm3: {big.sum()} (sum {r[big].sum():+.1f}); all {r.sum():+.1f}")
    for k in ('merge', 'detach', 'rift', 'subduction initiation'):
        m = np.array([x[1][k] > 0 for x in rows])
        print(f"   steps with {k}: {m.sum()} ({m.mean():.3f} of steps); of the big ones {m[big].mean():.2f}; residual there {r[m].sum():+.1f}, elsewhere {r[~m].sum():+.1f}")
