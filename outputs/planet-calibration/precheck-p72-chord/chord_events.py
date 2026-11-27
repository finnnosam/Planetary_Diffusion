"""Pass 72 (journal 95): the splits a run made, read from its stats.json event_list: rift events (count; new-plate area
p25 / p50 / p75; the share under 5 Mkm^2), how many were drawn as a 'chord' (rift_line 'chord': events carry 'line'),
the candidate lines' chord / bisector counts (stats chord_stats), merges and detachments, plates at 0 Myr. Pass 73
(journal 96): rift_line 'force' tags its cuts in tension 'force' (counted beside 'chord'); means print with one decimal.
Usage: chord_events.py <arm folder> [<arm folder> ...] (each holding seed*/stats.json) or seed folders."""
import json
import sys
from pathlib import Path

import numpy as np


def seeds(p):
    p = Path(p)
    return [p] if (p/'stats.json').exists() else sorted(q for q in p.glob('seed*') if (q/'stats.json').exists())


def read(sd):
    st = json.loads((sd/'stats.json').read_text())
    ev = st.get('event_list', [])
    rifts = [e for e in ev if e['kind'] == 'rift']
    a = np.array([e.get('area_mkm2', np.nan) for e in rifts], float)
    lines = [e.get('line') for e in rifts]
    pres = st.get('present', {})
    return {'rifts': len(rifts), 'area': a, 'chord': lines.count('chord'), 'bisector': lines.count('bisector'),
            'force': lines.count('force'),
            'merge': sum(e['kind'] == 'merge' for e in ev), 'detach': sum(e['kind'] == 'detach' for e in ev),
            'cand': st.get('chord_stats'), 'plates': pres.get('plates'), 'plates1': pres.get('plates_over_1mkm2')}


def line(r):
    a = r['area'][np.isfinite(r['area'])]
    q = np.percentile(a, [25, 50, 75]) if a.size else [np.nan]*3
    s = (f"rifts {r['rifts']:5.1f} (force {r['force']:.1f}, chord {r['chord']:.1f}, bisector {r['bisector']:.1f}) | new plate p25/50/75 "
         f"{q[0]:.1f} / {q[1]:.1f} / {q[2]:.1f} Mkm2, < 5 Mkm2 {np.mean(a < 5) if a.size else np.nan:.2f} | merges "
         f"{r['merge']:.1f} detach {r['detach']:.1f} | plates at 0 Myr {r['plates']:.1f} ({r['plates1']:.1f} >= 1 Mkm2)")
    if r['cand']:
        s += f" | candidates force {r['cand'].get('force', 0)} chord {r['cand'].get('chord')} bisector {r['cand'].get('bisector')}"
    return s


for arm in sys.argv[1:]:
    rs = []
    for sd in seeds(arm):
        r = read(sd); rs.append(r)
        print(f"{sd.parent.name}/{sd.name}: {line(r)}")
    if len(rs) > 1:
        m = {'rifts': np.mean([r['rifts'] for r in rs]), 'area': np.concatenate([r['area'] for r in rs]),
             'chord': np.mean([r['chord'] for r in rs]), 'bisector': np.mean([r['bisector'] for r in rs]),
             'force': np.mean([r['force'] for r in rs]),
             'merge': np.mean([r['merge'] for r in rs]), 'detach': np.mean([r['detach'] for r in rs]), 'cand': None,
             'plates': np.mean([r['plates'] or 0 for r in rs]), 'plates1': np.mean([r['plates1'] or 0 for r in rs])}
        print(f"{Path(arm).name} MEAN ({len(rs)} seeds; areas pooled): {line(m)}")
