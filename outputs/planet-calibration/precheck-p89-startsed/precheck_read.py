"""Pass 89 (journal 115): read the full-length cloud pre-check of the Start Sediment Sweep (seed 1018, unused; control
= today's PARAMS, 'all' = start_sediment 'all'). The start's sediment is a start condition, not a rule that accumulates
state; the read follows where the start's cover goes: the continents' sediment and its place against time, the land
and the sea, the continents on the datum (the spread's low side), the volume (crust + sediment, and the crust alone),
the highest node, plates, and the gates. Read only, not pooled.
Run from the repo root: precheck_read.py <runs folder> [seed] [arm ...]"""
import sys, glob, json, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_reviews import elevation_check as ec
from planet_reviews.crust_review import run_params
from planet_landscape import graph_landscape as gl
from pathlib import Path
g = CubeSphere(217); a = g.area/g.area.sum(); AREA = 4*np.pi*6371e3**2
seed = sys.argv[2] if len(sys.argv) > 2 else '1018'
arms = sys.argv[3:] or ['control', 'all']
def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); c = (c-.5*w[o])/c[-1]; return np.interp(q, c, x[o])
YR = 250e6
for name in arms:
    sd = f'{sys.argv[1]}/{name}/seed{seed}'
    if not glob.glob(sd+'/stats.json'):
        print(f'{name}: not finished'); continue
    st = json.load(open(f'{sd}/stats.json')); p = run_params(Path(sd)); ser = st['series']
    T = np.array([x['time_myr'] for x in ser]); pick = lambda key, t: ser[int(np.argmin(abs(T-t)))][key]
    print(f"{name}: plates at -200/-150/-100/-50/0 Myr " + ' / '.join(str(pick('plates', t)) for t in (-200, -150, -100, -50, 0))
          + f"; detached area max {max(x['detached_area_frac'] for x in ser):.4f}; ocean age mean at 0 {pick('ocean_age_mean_myr', 0):.0f} Myr")
    files = sorted(glob.glob(sd+'/snapshots/t*.npz'), key=lambda f: float(Path(f).stem[1:]))
    zmax = max((float(np.load(f)['elevation_m'].max()), Path(f).stem) for f in files)
    print(f"{name}: highest node over the saves {zmax[0]:.0f} m at {zmax[1]}")
    print('   t   sea  land  p25/p50/p90      room  cont.sed Mkm3 above/under sea  floor<500km/far  cont>=30km on datum p25/p50/p75  under  sed on low 0-20% of cont (bulk m)')
    for f in files[::max(1, len(files)//12)]+[files[-1]]:
        q = np.load(f); z = q['elevation_m'].astype(float); sl = float(q['sea_level_m']); zd = z+sl
        land = z > 0; s = q['sediment_bulk_m'].astype(float); Tc = q['crust_thickness_m'].astype(float)
        cont = Tc+s >= 20e3; c30 = Tc+s >= 30e3
        dist = gl.distance_to(cont, g.nbr, g.dist, 3e6)/1e3
        V = lambda m: (s[m]*a[m]).sum()*AREA/1e15
        cz = wq(zd[c30], a[c30], [.25, .5, .75])
        low = c30 & (zd <= wq(zd[c30], a[c30], .2))
        print(f"{Path(f).stem:>8} {sl:5.0f} {a[land].sum():.3f} " + '/'.join(f'{v:.0f}' for v in wq(z[land], a[land], [.25, .5, .9]))
              + f"  {np.sum(a*np.maximum(-zd, 0)):6.0f}  {V(cont & land):5.0f} / {V(cont & ~land):5.0f}        {V(~cont & (dist < 500)):5.0f} / {V(~cont & (dist >= 500)):4.0f}"
              f"       {cz[0]:5.0f} / {cz[1]:4.0f} / {cz[2]:4.0f}   {a[c30 & (zd < 0)].sum()/a[c30].sum():.3f}  {(s[low]*a[low]).sum()/a[low].sum():5.0f}")
    E = st['present']['erosion_totals_km3']; e = lambda k: E.get(k, 0.)/YR
    t = ec.snapshot_tests(files[-1], g, p)
    print(f"{name}: export {-e('continental_eroded_m3')-e('continental_deposited_m3'):.2f} km3/yr; volume change "
          f"{st['continental_volume_change']:+.3f} (crust {st['continental_crust_volume_change']:+.3f}); gates at 0 Myr: "
          + ', '.join(f'{k} {v:.3f}' if abs(v) < 5 else f'{k} {v:.0f}' for k, v in t.items() if isinstance(v, float)))
    print()
