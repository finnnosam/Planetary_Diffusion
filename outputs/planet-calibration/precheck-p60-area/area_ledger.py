# Pass 60: the planet's continental area ledger by process (flag area, stats.json series), Mkm2/Myr
import json, sys, numpy as np
S = 'outputs/planet-calibration/sweep-20261101-191306'
A = 510.0645  # Mkm2
KEYS = ('collision', 'subduction', 'rift', 'arc', 'holes', 'merge')
for combo in ('combo-02', 'combo-04'):
    rows = {}
    for seed in range(1001, 1007):
        d = json.load(open(f'{S}/{combo}/seed{seed}/stats.json'))
        s = d['series']; t = np.array([x['time_myr'] for x in s])
        L = {k: np.array([x['area_ledger_frac'].get(k, 0.) for x in s]) for k in KEYS}
        cf = np.array([x['continental_fraction'] for x in s])
        ch = np.array([x['continental_area_change_frac'] for x in s])
        i40 = int(np.argmin(abs(t+40)))
        for lab, i0 in (('250', 0), ('last40', i40)):
            dur = t[-1]-t[i0]
            r = {k: (L[k][-1]-L[k][i0])*A/dur for k in KEYS}
            r['sum'] = sum(r.values())
            r['flag_change'] = (cf[-1]-cf[i0])*A/dur
            r['area_change_frac'] = (ch[-1]-ch[i0])*A/dur
            r['other'] = r['flag_change']-r['sum']
            rows.setdefault(lab, []).append(r)
        rows.setdefault('cf', []).append((cf[0], cf[-1]))
    print(combo, 'flag continental fraction start/end (mean of six):', np.round(np.mean(rows['cf'], 0), 4))
    for lab in ('250', 'last40'):
        ks = list(rows[lab][0])
        print(f'  {lab:7s}' + ''.join(f'{k:>12s}' for k in ks))
        print('  mean   ' + ''.join(f'{np.mean([r[k] for r in rows[lab]]):12.3f}' for k in ks))
        print('  min    ' + ''.join(f'{np.min([r[k] for r in rows[lab]]):12.3f}' for k in ks))
        print('  max    ' + ''.join(f'{np.max([r[k] for r in rows[lab]]):12.3f}' for k in ks))
