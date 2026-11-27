"""Pass 87 (journal 113): read the full-length cloud pre-check of the Arc Accretion Sweep's docking arms (seed 1017,
unused; today's PARAMS = control; bd300 = trench_polarity 'pair' + trench_flux 'consumed' + arc_accretion, arc_H_ocean_m
300; bd1000 = the same with arc_H_ocean_m 1,000). Docking accumulates state and ran away on one seed on 27 Oct, so
the read is first the runaway (plates, detached area, trench against time), then the arcs' budget, the land, the sea,
the crust, the highest node and the gates at 0 Myr. Read only, not pooled.
Run from the repo root: precheck_read.py <runs folder> [seed] [arm ...]"""
import sys, glob, json, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_reviews import elevation_check as ec
g = CubeSphere(217); a = g.area/g.area.sum(); AREA = 4*np.pi*6371e3**2
seed = sys.argv[2] if len(sys.argv) > 2 else '1017'
arms = sys.argv[3:] or ['control', 'bd300', 'bd1000']
def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); c = (c-.5*w[o])/c[-1]; return np.interp(q, c, x[o])
YR = 250e6
for name in arms:
    sd = f'{sys.argv[1]}/{name}/seed{seed}'
    if not glob.glob(sd+'/stats.json'):
        print(f'{name}: not finished'); continue
    st = json.load(open(f'{sd}/stats.json'))
    p = {**P.PARAMS, **{k: v for k, v in st['params_effective'].items() if k in P.PARAMS}}
    ser = st['series']
    T = np.array([x['time_myr'] for x in ser])
    pick = lambda key, t: ser[int(np.argmin(abs(T-t)))][key]
    print(f"{name}: plates at -200/-150/-100/-50/0 Myr " + ' / '.join(str(pick('plates', t)) for t in (-200, -150, -100, -50, 0))
          + f"; max {max(x['plates'] for x in ser)}; detached area max {max(x['detached_area_frac'] for x in ser):.4f}; "
          f"trench (k km) at -200/-150/-100/-50/0 " + ' / '.join(f"{pick('subduction_km', t)/1e3:.1f}" for t in (-200, -150, -100, -50, 0))
          + f"; ridge at 0 {pick('ridge_km', 0)/1e3:.1f}k; ocean age mean at 0 {pick('ocean_age_mean_myr', 0):.0f} Myr; "
          f"mean speed -100..0 {np.mean([np.mean(x['speeds_cm_yr']) for x in ser if x['time_myr'] >= -100]):.2f} cm/yr")
    ev = st['event_list']
    print(f"{name}: events rift (continental) {sum(1 for e in ev if e.get('kind') == 'rift' and 'line' in e)}, "
          f"rift (ocean) {sum(1 for e in ev if e.get('kind') == 'rift' and 'line' not in e)}, "
          f"initiation {sum(1 for e in ev if e.get('kind') == 'subduction initiation')}, detach {sum(1 for e in ev if e.get('kind') == 'detach')}, "
          f"merge {sum(1 for e in ev if e.get('kind') == 'merge')}")
    files = sorted(glob.glob(sd+'/snapshots/t*.npz'))
    zmax = max((float(np.load(f)['elevation_m'].max()), f[-12:-4]) for f in files)
    for tag in ('t-0100.0', 't-0040.0', 't+0000.0'):
        q = np.load(f'{sd}/snapshots/{tag}.npz'); z = q['elevation_m'].astype(float); sl = float(q['sea_level_m']); zd = z+sl
        land = z > 0; s = q['sediment_bulk_m'].astype(float); Tc = q['crust_thickness_m'].astype(float)+s
        cont = Tc >= 20e3
        ca = a[cont].sum()
        print(f"{name} {tag}: sea {sl:.0f} m; land {a[land].sum():.3f}, p25/p50/p90 " + '/'.join(f'{v:.0f}' for v in wq(z[land], a[land], [.25, .5, .9]))
              + f"; room {np.sum(a*np.maximum(-zd, 0)):.0f} m; continent {ca:.3f}, 20-30 km {a[cont & (Tc < 30e3)].sum():.3f} of the surface, "
              f"> 45 / > 55 km {a[cont & (Tc > 45e3)].sum()/ca:.3f} / {a[cont & (Tc > 55e3)].sum()/ca:.3f}, crust p50 {wq(Tc[cont]/1e3, a[cont], .5):.1f} km")
    E = st['present']['erosion_totals_km3']
    e = lambda k: E.get(k, 0.)/YR
    t = ec.snapshot_tests(f'{sd}/snapshots/t+0000.0.npz', g, p)
    print(f"{name}: highest node over the saves {zmax[0]:.0f} m at {zmax[1]}; arcs (km3/yr, whole run): growth on ocean arcs "
          f"{e('subduction_on_ocean_arcs_m3'):.2f}, accreted {e('arc_accreted_m3'):.2f} (on continents {e('arc_accreted_on_continent_m3'):.2f}), "
          f"docked and spread on oceanic plates {e('arc_docked_spread_m3'):.2f}; export {-e('continental_eroded_m3')-e('continental_deposited_m3'):.2f}; "
          f"volume change {st['continental_volume_change']:+.3f} (crust {st['continental_crust_volume_change']:+.3f}); gates at 0 Myr: "
          + ', '.join(f'{k} {v:.3f}' if abs(v) < 5 else f'{k} {v:.0f}' for k, v in t.items() if isinstance(v, float)))
    print()
