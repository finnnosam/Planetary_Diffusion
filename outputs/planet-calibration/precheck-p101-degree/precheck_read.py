"""Pass 99 pre-check reader (cloud, seed 1024, planet only; adapted from pass 97's): per arm folder, the gates pooled over the last 40 Myr,
the highest node in the saved snapshots, land share and p25 / p50 / p90 as run, the sea on the datum, the
continents (crust + bulk sediment >= 30 km... crystalline >= 30 km as continent.py) on the datum, the drowned band and
the submerged continent's depth by distance from land (shelf_profile's bins) at 0 Myr."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'precheck-p97-tail'))
from drowned_band import rows, BANDS                               # noqa: E402
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import run_params, sediment_bulk_m
from planet_reviews.elevation_check import snapshot_tests, weighted_median
g = CubeSphere(217); a = g.area
BINS = (0, 50, 100, 150, 200, 300, 400, 600, 1000)


def wq(v, w, q):
    o = np.argsort(v); c = np.cumsum(w[o])
    return float(v[o][np.searchsorted(c, q*c[-1])]) if c.size and c[-1] > 0 else np.nan


for arm in map(Path, sys.argv[1:]):
    sd = next(arm.glob('seed*')); p = run_params(sd)
    snaps = sorted((sd/'snapshots').glob('t*.npz'), key=lambda s: float(s.stem[1:]))
    last = [s for s in snaps if float(s.stem[1:]) >= -40.5]
    tests = [snapshot_tests(s, g, p) for s in last]
    o = {k: weighted_median(np.concatenate([t['_'+k][0] for t in tests]), np.concatenate([t['_'+k][1] for t in tests]))
         for k in ('young_ocean', 'thin_crust')}
    zmax = max((float(np.load(s)['elevation_m'].max()), s.stem) for s in snaps)
    print(f"== {arm.name}/{sd.name}: gates pooled young-ocean {o['young_ocean']:+.0f} m, thin-crust {o['thin_crust']:+.0f} m; "
          f"highest node in the snapshots {zmax[0]:.0f} m ({zmax[1]})")
    for t in ('t-0250.0', 't+0000.0'):
        f = np.load(sd/'snapshots'/f'{t}.npz')
        Tc = f['crust_thickness_m'].astype(float); Sb = sediment_bulk_m(f); z = f['elevation_m'].astype(float)
        sea = float(f['sea_level_m']); zd = z+sea
        land = z > 0
        c30 = Tc >= 30e3
        qs = [wq(zd[c30], a[c30], q) for q in (.1, .25, .5, .75, .9)]
        print(f"  {t}: land {a[land].sum()/a.sum():.3f}, p25 / p50 / p90 {wq(z[land], a[land], .25):.0f} / "
              f"{wq(z[land], a[land], .5):.0f} / {wq(z[land], a[land], .9):.0f} m; sea on the datum {sea:+.0f} m; continents "
              f"(crust >= 30 km) on the datum p10-p90 " + ' / '.join(f'{x:.0f}' for x in qs) +
              f", under the datum {a[c30 & (zd < 0)].sum()/a[c30].sum():.3f}")

        dyn = f['dynamic_topography_m'].astype(float); air = f['dynamic_topography_air_m'].astype(float)
        cc = Tc >= 30e3
        q6 = (.05, .1, .25, .5, .75, .9)
        print(f"    field: rms air-loaded {np.sqrt((a*air**2).sum()/a.sum()):.0f} m; on continents (crystalline >= 30 km) "
              f"under / above the datum {(a*dyn)[cc & (zd < 0)].sum()/a[cc & (zd < 0)].sum():+.0f} / "
              f"{(a*dyn)[cc & (zd >= 0)].sum()/a[cc & (zd >= 0)].sum():+.0f} m; land {(a*dyn)[land].sum()/a[land].sum():+.0f}, "
              f"ocean {(a*dyn)[~land].sum()/a[~land].sum():+.0f} m")
        print('    continents on the datum z p5..p90: ' + ' / '.join(f'{wq(zd[cc], a[cc], q):.0f}' for q in q6) +
              '; z - dyn: ' + ' / '.join(f'{wq((zd-dyn)[cc], a[cc], q):.0f}' for q in q6))
        if 'slab_mass_young_kg' in f.files:
            from planet_landscape import graph_landscape as gl2
            sm = f['slab_mass_young_kg'].astype(float)
            ds = gl2.distance_to(sm > 0, g.nbr, g.dist)/1e3
            print(f"    corr(air field, distance to young slab mass) {np.corrcoef(air, ds)[0, 1]:+.2f}; "
                  f"raw flow power l1-30 rms {np.sqrt(f['flow_power_raw_m2'][1:31].sum()/(4*np.pi)):.0f} m")
        age = f['ocean_age_myr'].astype(float)
        print(f"    young floor 0-60 Myr {a[np.isfinite(age) & (age < 60)].sum()/a.sum():.3f} of the surface; ocean age "
              f"mean {(a*np.nan_to_num(age))[np.isfinite(age)].sum()/a[np.isfinite(age)].sum():.0f} Myr")
        rs = rows(g, Tc, Sb, z)
        print('    continent by elevation about the sea: ' + ' | '.join(
            f"{lo:.0f}..{hi:.0f} {r['share']:.3f} (Tc {r['Tc']:.1f} km, S {r['S']:.0f} m, {r['dist']:.0f} km in)"
            for (lo, hi), r in zip(BANDS, rs)))
        if t == 't+0000.0':
            cont = (Tc+Sb) >= 20e3; sub = cont & (z < 0)
            d = gl.distance_to(z >= 0, g.nbr, g.dist)/1e3
            out = []
            for lo, hi in zip(BINS[:-1], BINS[1:]):
                m = sub & (d >= lo) & (d < hi)
                out.append(f'{lo}-{hi}: {a[m].sum()/a[sub].sum():.2f} {wq(z[m], a[m], .25):.0f}/{wq(z[m], a[m], .5):.0f}/'
                           f'{wq(z[m], a[m], .75):.0f}')
            print(f'    submerged continent {a[sub].sum()/a.sum():.3f} by km from land (share, p25/p50/p75): ' + ' | '.join(out))
