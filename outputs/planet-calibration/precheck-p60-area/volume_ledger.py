# Pass 60: the planet's continental volume terms over 250 Myr (erosion_totals and the crust/sediment ledgers), km3/yr
import json, numpy as np
S = 'outputs/planet-calibration/sweep-20261101-191306'
for combo in ('combo-02', 'combo-04'):
    E = []; L = []; V = []
    for seed in range(1001, 1007):
        d = json.load(open(f'{S}/{combo}/seed{seed}/stats.json')); s = d['series']
        last = d['present']; assert last['time_myr'] == 0.0; E.append(last['erosion_totals_km3']); L.append(last['sediment_ledger_km3'])
        V.append((s[0]['continental_volume_km3'], last['continental_volume_km3'], s[0]['continental_sediment_km3'], last['continental_sediment_km3'], d['continental_volume_change'], d['continental_crust_volume_change']))
    yrs = 250e6
    print(combo, 'erosion_totals (km3/yr over 250 Myr; its _m3 keys hold km3, six-seed mean):')
    for k in E[0]:
        v = np.mean([e.get(k, 0.) for e in E])
        if k.endswith('_m3'): print('   %-40s %8.3f' % (k, v/250e6))
    V = np.array(V); print('   continental volume start / end (Mkm3) %.0f / %.0f; sediment %.0f / %.0f; change total %.3f crust %.3f' % (*(V[:, :4].mean(0)/1e6), V[:, 4].mean(), V[:, 5].mean()))
