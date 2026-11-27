"""Pass 50, reported only (P50_RULE): where on Earth T1 (E3 topography 0.15-0.30 km) and T2 (geoid VR l 2-8 >= 0.5,
both corrections) hold over a two-parameter family: lithosphere 0-100 km 10x, asthenosphere 100-300 km at a contrast,
upper mantle 300-670 km 1, lower mantle 670-2891 km at a ratio. A measurement of the trade-off, not a candidate set."""
import json
import numpy as np
import profiles as p

rows = []
print('LM ratio | asth: topo m / min VR (l2-8) / corr')
for lm in (3., 10., 30., 100., 300.):
    line = []
    for ast in (1., .3, .1, .03, .01):
        base, eta = p.layered([(2891, lm), (670, 1.), (300, ast), (100, 10.)])
        r = p.evaluate(base, eta, full=False)[0]
        v = min(r['geoid'][k]['VR_l2_8'] for k in p.OBS); c = min(r['geoid'][k]['corr_l2_8'] for k in p.OBS)
        ok = (150 <= r['topo_rms_m'] <= 300) and v >= .5
        rows.append({'lm': lm, 'asth': ast, 'topo_rms_m': r['topo_rms_m'], 'min_VR': v, 'min_corr': c,
                     'amp': r['geoid']['chambat']['amp_ratio_l2_8'], 'both': ok})
        line.append(f"{r['topo_rms_m']:4.0f}/{v:+.2f}/{c:+.2f}{'*' if ok else ' '}")
    print(f'{lm:5.0f}  | ' + ' | '.join(f'a{a:g}: {s}' for a, s in zip((1, .3, .1, .03, .01), line)))
json.dump(rows, open('scan.json', 'w'), indent=1)
