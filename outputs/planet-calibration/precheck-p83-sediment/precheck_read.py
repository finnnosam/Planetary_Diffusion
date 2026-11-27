"""Pass 83 (journal 108): read the full-length cloud pre-check of S (seed 1014, unused; today's PARAMS, the control and
margin_wedge True). Land share and p25 / p50 / p90, sea on the datum, room below the datum, the floor's sediment by
distance from the continent (sediment_where's classes), the highest node over the saves, the elevation gates at 0 Myr,
the wedge's booked volume. Read only, not pooled. Run from the repo root: precheck_read.py <runs folder>"""
import sys, glob, json, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import plates as P
from planet_reviews import elevation_check as ec
g = CubeSphere(217); a = g.area/g.area.sum(); AREA = 4*np.pi*6371e3**2
def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); c = (c-.5*w[o])/c[-1]; return np.interp(q, c, x[o])
for name, extra in (('control', {}), ('wedge', {'margin_wedge': True})):
    sd = f'{sys.argv[1]}/{name}/seed1014'; p = {**P.PARAMS, **extra}
    files = sorted(glob.glob(sd+'/snapshots/t*.npz'))
    zmax = max((float(np.load(f)['elevation_m'].max()), f[-12:-4]) for f in files)
    for tag in ('t-0250.0', 't-0100.0', 't-0040.0', 't+0000.0'):
        q = np.load(f'{sd}/snapshots/{tag}.npz'); z = q['elevation_m'].astype(float); sl = float(q['sea_level_m']); zd = z+sl
        land = z > 0; s = q['sediment_bulk_m'].astype(float); T = q['crust_thickness_m'].astype(float)+s
        cont = T >= 20e3; dist = gl.distance_to(cont, g.nbr, g.dist, 3e6)/1e3
        V = lambda m: float((s*a)[m].sum()*AREA/1e15)
        print(f"{name} {tag}: sea {sl:.0f} m; land {a[land].sum():.3f}, p25/p50/p90 " + '/'.join(f'{v:.0f}' for v in wq(z[land], a[land], [.25, .5, .9]))
              + f"; room {np.sum(a*np.maximum(-zd, 0)):.0f} m; sediment Mkm3: continent above sea {V(cont & land):.1f}, under sea {V(cont & ~land):.1f}, "
              f"floor 0-200 {V(~cont & (dist < 200)):.1f}, 200-500 {V(~cont & (dist >= 200) & (dist < 500)):.1f}, 500-1000 {V(~cont & (dist >= 500) & (dist < 1000)):.1f}, > 1000 {V(~cont & (dist >= 1000)):.1f}")
    st = json.load(open(f'{sd}/stats.json')); E = st['present']['erosion_totals_km3']
    t = ec.snapshot_tests(f'{sd}/snapshots/t+0000.0.npz', g, p)
    print(f"{name}: highest node over the saves {zmax[0]:.0f} m at {zmax[1]}; wedge booked {E.get('margin_wedge_m3', 0.)/250e6:.3f} km3/yr "
          f"(spill path {E.get('margin_wedge_spill_m3', 0.)/250e6:.3f}), deep fill {E.get('deep_sea_sediment_m3', 0.)/250e6:.3f}, to the floor "
          f"{E.get('ocean_sediment_m3', 0.)/250e6:.3f}, subducted {E.get('sediment_subducted_m3', 0.)/250e6:.3f}; volume change {st['continental_volume_change']:+.3f}; gates at 0 Myr: "
          + ', '.join(f'{k} {v:.3f}' if abs(v) < 5 else f'{k} {v:.0f}' for k, v in t.items() if isinstance(v, float)))
