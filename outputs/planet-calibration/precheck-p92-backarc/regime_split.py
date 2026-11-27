"""Pass 92: the upper plate's regime on continental overriding trench (u - rollback: shortening > 0, extension < 0;
output-only record 'overriding_*' in every run since pass 92) and the back-arc steps, per arm and seed, whole run and
the last 100 Myr, against Earth's transects (Lallemand et al. 2008, Table 1, 166 transects at 2 deg; continental upper
plates by name, 126: extension 0.11, neutral 0.67, compression 0.21; all: 0.27 / 0.52 / 0.22; journal 119 sec. 2).
Args: arm folders (each with seed*/stats.json) or run folders. Run from anywhere."""
import json, sys
from pathlib import Path
import numpy as np


def window(series, t0):
    a = next(r for r in series if r['time_myr'] >= t0-1e-6)['erosion_totals_km3']
    b = series[-1]['erosion_totals_km3']
    return {k: b.get(k, 0.)-a.get(k, 0.) for k in b}


for arm in sys.argv[1:]:
    rows = []
    for f in sorted(Path(arm).glob('seed*/stats.json')):
        S = json.loads(f.read_text())
        ser = S['series']
        r = [S['seed']]
        for t0 in (ser[0]['time_myr'], -100.):
            E = window(ser, t0)
            L = max(E.get('overriding_cont_trench_m', 0.), 1e-30)
            yrs = (ser[-1]['time_myr']-t0)*1e6
            r += [E.get('overriding_cont_trench_m', 0.)/max(E.get('overriding_trench_m', 0.), 1e-30),
                  E.get('overriding_cont_shortening_m', 0.)/L, E.get('overriding_cont_extending_m', 0.)/L,
                  E.get('overriding_cont_shortening5_m', 0.)/L, E.get('overriding_cont_extending5_m', 0.)/L,
                  E.get('backarc_steps', 0.)*1e9/(yrs/1e6), E.get('backarc_stretched_m3', 0.)/yrs,
                  E.get('backarc_area_m2', 0.)/max(E.get('backarc_expected_m2', 0.), 1e-30) if E.get('backarc_expected_m2') else np.nan]
        rows.append(r)
    if not rows:
        continue
    R = np.array(rows, float)
    print(f'\n{Path(arm).name}: continental overriding trench, share of all overriding / shortening / extending / '
          f'shortening > 5 mm/yr / extending > 5 mm/yr; back-arc steps per Myr, stretched km3/yr, area taken / expected')
    for lab, o in (('whole run', 1), ('last 100 Myr', 9)):
        for row in R:
            print(f'  {lab:12s} seed {int(row[0])}: ' + ' / '.join(f'{x:.2f}' for x in row[o:o+5]) +
                  f';  {row[o+5]:.1f} / {row[o+6]:.3f} / {row[o+7]:.2f}')
        m = np.nanmean(R[:, o:o+8], 0)
        print(f'  {lab:12s} mean     : ' + ' / '.join(f'{x:.2f}' for x in m[:5]) + f';  {m[5]:.1f} / {m[6]:.3f} / {m[7]:.2f}')
print('\nEarth (Lallemand 2008, continental upper plates): extension 0.11, neutral 0.67, compression 0.21 of transects')
