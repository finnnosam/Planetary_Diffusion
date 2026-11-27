"""Pass 103: read split_log_run.py's log. Per 50 Myr window: for the plate holding the most continent each step (the
supercontinent plate) and for all other plates with continent: steps, the line's kind (force / chord / bisector), the
median L and part area, the mean acceptance probability, events made, the margin and continental area shed, and the
boundary make-up of the big plate (share of its edge nodes trench / overriding / ridge)."""
import json, sys
from collections import defaultdict
import numpy as np
for path in sys.argv[1:]:
    log = json.load(open(path))
    byt = defaultdict(list)
    for r in log: byt[r['t']].append(r)
    print(path)
    for t0 in range(-250, 0, 50):
        big, rest = [], []
        for t, rs in byt.items():
            if t0 <= t < t0+50:
                b = max(rs, key=lambda r: r['Ac']); big.append(b); rest += [r for r in rs if r is not b]
        for name, rs in (('big', big), ('others', rest)):
            if not rs: continue
            k = defaultdict(int)
            for r in rs: k[r['kind']] += 1
            made = [r for r in rs if r['made']]
            line = (f'  {t0:4d}..{t0+50:<4d} {name:6} n {len(rs):4d} Ac med {np.median([r["Ac"] for r in rs]):6.1f} | kinds ' +
                    ' '.join(f'{kk} {v/len(rs):.2f}' for kk, v in sorted(k.items())) +
                    f' | L med {np.median([r["L"] for r in rs]):6.0f} km, part med {np.median([r["part"] for r in rs]):5.1f} (cont {np.median([r["part_c"] for r in rs]):5.1f}) | p mean {np.mean([r["p"] for r in rs]):.3f}'
                    f' | made {len(made)}: L sum {sum(r["L"] for r in made):6.0f} km, cont shed {sum(r["part_c"] for r in made):5.1f}')
            if name == 'big':
                e = np.array([[r['edge'], r['trench'], r['over'], r['ridge']] for r in rs], float)
                line += f' | edge trench / over / ridge {e[:,1].sum()/e[:,0].sum():.2f} / {e[:,2].sum()/e[:,0].sum():.2f} / {e[:,3].sum()/e[:,0].sum():.2f}'
                for kk in ('force', 'chord', 'bisector'):
                    s = [r for r in rs if r['kind'] == kk]
                    if s: line += f'\n               {kk:8}: L med {np.median([r["L"] for r in s]):6.0f}, part_c med {np.median([r["part_c"] for r in s]):5.1f}, p {np.mean([r["p"] for r in s]):.3f}, made {sum(r["made"] for r in s)}'
            print(line)
