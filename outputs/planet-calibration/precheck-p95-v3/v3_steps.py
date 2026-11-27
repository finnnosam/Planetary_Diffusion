"""Per-step crust residual (km3) from a run's stats series: crust volume change less the booked terms (pass 95).
Prints the residual's total, the share in the largest 10 steps, and those steps. Terms as crust_budget.budget
(the arc conversion term from the area ledger)."""
import json, sys
from pathlib import Path
import numpy as np
from planet_reviews.crust_review import run_params
R = 6371.0
KEYS = ['continental_crust_eroded_m3', 'subduction_on_continents_m3', 'sediment_accreted_m3',
        'collision_sediment_to_crust_m3', 'volcanic_added_m3', 'arc_accreted_on_continent_m3']
for seed in sys.argv[1:]:
    st = json.loads((Path(seed)/'stats.json').read_text()); p = run_params(seed)
    A = 4*np.pi*R**2
    S = st['series']; res = []
    for a, b in zip(S[:-1], S[1:]):
        E = lambda k: b['erosion_totals_km3'].get(k, 0.)-a['erosion_totals_km3'].get(k, 0.)
        arc = (b['area_ledger_frac'].get('arc', 0.)-a['area_ledger_frac'].get('arc', 0.))*A*p['arc_to_continent_thick_m']/1e3
        terms = sum(E(k) for k in KEYS)+arc-E('eclogite_foundered_m3')
        res.append((b['time_myr'], b['continental_volume_km3']-a['continental_volume_km3']-terms))
    r = np.array([x[1] for x in res]); t = np.array([x[0] for x in res])
    o = np.argsort(-np.abs(r))[:10]
    last = t > -40
    print(f"{Path(seed).parent.name}/{Path(seed).name}: total {r.sum()/1e6:+.1f} Mkm3 (last 40 {r[last].sum()/1e6:+.1f}); "
          f"largest 10 steps {r[o].sum()/1e6:+.1f}; rest {(r.sum()-r[o].sum())/1e6:+.1f}; median |step| {np.median(np.abs(r))/1e6:.2f}")
    print('   ', ', '.join(f"{t[i]:.0f}: {r[i]/1e6:+.1f}" for i in o))
