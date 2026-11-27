"""Pass 85 (journal 110): read the full-length cloud pre-check of sediment_entrainment (seed 1016, unused; today's
PARAMS, crust_flow_rule 'rey', with sediment_entrainment bedrock / cover). Land share and p25 / p50 / p90, sea on the datum, room
below the datum, label-free continent and its > 45 / > 55 km shares, sediment above the sea and under it, the highest
node over the saves, the elevation gates at 0 Myr, export and erosion, the flow's substeps. Read only, not pooled.
Run from the repo root: precheck_read.py <runs folder> [seed]"""
import sys, glob, json, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_reviews import elevation_check as ec
g = CubeSphere(217); a = g.area/g.area.sum(); AREA = 4*np.pi*6371e3**2
seed = sys.argv[2] if len(sys.argv) > 2 else '1016'
def wq(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); c = (c-.5*w[o])/c[-1]; return np.interp(q, c, x[o])
for name in ('bedrock', 'cover'):
    sd = f'{sys.argv[1]}/{name}/seed{seed}'
    if not glob.glob(sd+'/stats.json'):
        print(f'{name}: not finished'); continue
    p = {**P.PARAMS, 'sediment_entrainment': name}
    files = sorted(glob.glob(sd+'/snapshots/t*.npz'))
    zmax = max((float(np.load(f)['elevation_m'].max()), f[-12:-4]) for f in files)
    for tag in ('t-0200.0', 't-0100.0', 't-0040.0', 't+0000.0'):
        q = np.load(f'{sd}/snapshots/{tag}.npz'); z = q['elevation_m'].astype(float); sl = float(q['sea_level_m']); zd = z+sl
        land = z > 0; s = q['sediment_bulk_m'].astype(float); T = q['crust_thickness_m'].astype(float)+s
        cont = T >= 20e3
        V = lambda m: float((s*a)[m].sum()*AREA/1e15)
        ca = a[cont].sum()
        print(f"{name} {tag}: sea {sl:.0f} m; land {a[land].sum():.3f}, p25/p50/p90 " + '/'.join(f'{v:.0f}' for v in wq(z[land], a[land], [.25, .5, .9]))
              + f"; room {np.sum(a*np.maximum(-zd, 0)):.0f} m; continent {ca:.3f}, > 45 / > 55 km {a[cont & (T > 45e3)].sum()/ca:.3f} / "
              f"{a[cont & (T > 55e3)].sum()/ca:.3f}, crust p50 {wq(T[cont]/1e3, a[cont], .5):.1f} km; sediment Mkm3 continent above sea "
              f"{V(cont & land):.1f}, under sea {V(cont & ~land):.1f}, floor {V(~cont):.1f}")
    st = json.load(open(f'{sd}/stats.json')); E = st['present']['erosion_totals_km3']
    t = ec.snapshot_tests(f'{sd}/snapshots/t+0000.0.npz', g, p)
    yrs = 250e6
    print(f"{name}: highest node over the saves {zmax[0]:.0f} m at {zmax[1]}; eroded {E.get('continental_eroded_m3', 0.)/yrs:.2f} km3/yr, "
          f"deposited {E.get('continental_deposited_m3', 0.)/yrs:.2f} (submerged {E.get('continental_deposited_submerged_m3', 0.)/yrs:.2f}), crystalline {E.get('continental_crust_eroded_m3', 0.)/yrs:.2f}, "
          f"to the floor {E.get('ocean_sediment_m3', 0.)/yrs:.3f}, deep-sea fill {E.get('deep_sea_sediment_m3', 0.)/yrs:.3f}, foundered {-E.get('eclogite_foundered_m3', 0.)/yrs:.3f}; "
          f"flow substeps per Myr {E.get('crust_flow_substeps', 0.)*1e9/250:.0f}, edges at the ceiling "
          f"{E.get('crust_flow_capped_edges', 0.)/max(E.get('crust_flow_edges', 0.), 1e-30):.4f}, pinned share of the continent "
          f"{E.get('crust_flow_pinned_area', 0.)*1e9/250:.3f}; volume change {st['continental_volume_change']:+.3f}; gates at 0 Myr: "
          + ', '.join(f'{k} {v:.3f}' if abs(v) < 5 else f'{k} {v:.0f}' for k, v in t.items() if isinstance(v, float)))
