"""Pass 82 step-back, probe 9: the sea and the land through the run (stats.json series, six seeds), and each seed's
sediment budget over the run (km3/yr). Usage: python series.py [arm]"""
import json, sys, numpy as np
A = sys.argv[1] if len(sys.argv) > 1 else 'outputs/planet-calibration/sweep-20261109-222925/lith_frame-lamb'
S = {}
for s in range(1001, 1007):
    S[s] = {round(r['time_myr']): r for r in json.load(open(f'{A}/seed{s}/stats.json'))['series']}
print('time Myr | sea on the datum mean (sd) | land share | land p50 mean (sd) | land p25 | p90 | continental share | ocean age mean | continental sediment Mkm3')
for t in [t for t in sorted(S[1001]) if t % 25 == 0]:
    rows = [S[s][t] for s in S if t in S[s]]
    sea = np.array([r['sea_level_on_datum_m'] for r in rows]); p = np.array([[r['land_elevation_percentiles_m'][k] for k in ('25', '50', '90')] for r in rows])
    print(f"{t:5d} | {sea.mean():6.0f} ({sea.std(ddof=1):3.0f}) | {np.mean([r['land_fraction'] for r in rows]):.3f} | {p[:, 1].mean():5.0f} ({p[:, 1].std(ddof=1):3.0f}) | {p[:, 0].mean():5.0f} | "
          f"{p[:, 2].mean():5.0f} | {np.mean([r['continental_fraction'] for r in rows]):.3f} | {np.mean([r['ocean_age_mean_myr'] for r in rows]):5.1f} | {np.mean([r['continental_sediment_km3'] for r in rows])/1e6:6.1f}")
print('\nper seed, whole run (km3/yr): sea at 0 Myr | land p50 | eroded off continents | deposited on continents (of it under the sea) | terrigenous to the ocean floor | pelagic')
for s in S:
    j = json.load(open(f'{A}/seed{s}/stats.json')); e = j['present']['erosion_totals_km3']; yr = 250e6
    print(f"  {s} | {j['present']['sea_level_on_datum_m']:.0f} | {j['present']['land_elevation_percentiles_m']['50']:.0f} | {-e['continental_eroded_m3']/yr:.2f} | "
          f"{e['continental_deposited_m3']/yr:.2f} ({e['continental_deposited_submerged_m3']/yr:.2f}) | {e['ocean_sediment_m3']/yr:.2f} | {e['pelagic_m3']/yr:.2f}")
