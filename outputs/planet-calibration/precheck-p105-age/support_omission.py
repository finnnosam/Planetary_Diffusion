"""Pass 105 (journal 142): how much continental lithospheric support changes along the landscape's parcel paths, which the
landscape forcing does not carry (the review's item 2, docs/review/2026-11-review.md).

The landscape stage (planet_landscape.coupling.build_forcing) forces a continental parcel with the plate layer's tectonic
thickening (a rock load it compensates itself) and the change of dynamic topography less the sea's change. The plate
layer's elevation also has C, the mantle lithosphere's thermal support (lithosphere.elevation_term of the
tectonothermal age, `lithosphere` 'thermal_age'), which moves when the age changes: cooling (older, thicker, lower),
stretching (younger, thinner, higher), collision (thicker, lower), resets. The landscape never sees that change.

Parcels: the n = 217 nodes at 0 Myr (one per plate-layer node), walked back with coupling.step_back through the saved
every-Myr snapshots over the window. Read on parcels continental at every sample of the window (the forcing's own
continental channel), area-weighted:
  d_lith   C at 0 Myr minus C at the window's start, air-loaded (m); nearest-node lithosphere age (no blending of NaN)
  airy     sum of tectonic thickening x (1 - rho_crust / rho_mantle): the forcing's load at Airy equilibrium (m)
  dyn      sum of the forcing's surface channel on continent (dynamic topography change less the sea's change) (m)
  dz       the plate layer's own elevation change along the path (m); gap = dz - airy - dyn holds the plate layer's
           erosion and deposition, sediment loads, flexure, density and C together
Groups: all such parcels; land at 0 Myr; land above 3 km at 0 Myr; no thickening on the path (|sum| < 50 m: support
changing without crustal thickening, the review's test case); thickened (> 500 m); thinned (< -500 m).

  python support_omission.py <plate seed folder> [window_myr] > support_omission.txt
"""
import json
import sys
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import coupling
from planet_tectonics.lithosphere import elevation_term
from planet_tectonics.plates import PARAMS


def wq(v, w, qs=(10, 50, 90)):
    o = np.argsort(v); c = np.cumsum(w[o]); c /= c[-1]
    return [float(v[o][np.searchsorted(c, q/100)]) for q in qs]


def main(seed_dir, window=30):
    seed_dir = Path(seed_dir)
    stats = json.loads((seed_dir/'stats.json').read_text())
    p = {**PARAMS, **stats.get('params_effective', {})}
    if p.get('lithosphere') != 'thermal_age':
        print(f"lithosphere {p.get('lithosphere')!r}: no C term in this run"); return
    hist = coupling.PlateHistory(seed_dir/'snapshots')
    g = hist.g
    rc, rm = p['rho_crust'], p['rho_mantle']
    pos = g.xyz.copy()
    area = g.area.copy()

    def at(t, q):
        f, _ = hist.fields(float(t))
        i = fields.nearest_nodes(g, q)
        return f, i

    f0, i0 = at(0, pos)
    cont_all = f0['continental'][i0].astype(bool)
    lith0 = elevation_term(f0['lithosphere_age_myr'][i0].astype(np.float64), p)
    z0 = f0['elevation_m'][i0].astype(np.float64)
    airy = np.zeros(g.size); dyn = np.zeros(g.size)
    q = pos
    f_t, i_t = f0, i0
    for t in range(0, -window, -1):
        qp = coupling.step_back(hist, float(t), q)
        f_p, i_p = at(t-1, qp)
        cont_all &= f_p['continental'][i_p].astype(bool)
        airy += f_t['tectonic_thickening_m'][i_t]*(1-rc/rm)
        dyn += (f_t['dynamic_topography_m'][i_t]-f_p['dynamic_topography_m'][i_p]) - \
               (float(f_t['sea_level_m'])-float(f_p['sea_level_m']))
        q, f_t, i_t = qp, f_p, i_p
    lith_s = elevation_term(f_t['lithosphere_age_myr'][i_t].astype(np.float64), p)
    zs = f_t['elevation_m'][i_t].astype(np.float64)
    d_lith = lith0-lith_s
    dz = z0-zs
    gap = dz-airy-dyn
    thick_sum = airy/(1-rc/rm)
    groups = {'continental throughout': cont_all,
              'land at 0 Myr': cont_all & (z0 > 0),
              'land > 3 km at 0 Myr': cont_all & (z0 > 3000),
              'no thickening on the path (|sum| < 50 m), land': cont_all & (z0 > 0) & (np.abs(thick_sum) < 50),
              'thickened (> 500 m), land': cont_all & (z0 > 0) & (thick_sum > 500),
              'thinned (< -500 m), land': cont_all & (z0 > 0) & (thick_sum < -500)}
    print(f'{seed_dir}: window {window} Myr, parcels = n {g.n} nodes at 0 Myr; C = lithosphere.elevation_term (air-loaded)')
    print(f"lith_frame {p.get('lith_frame')}; rho_crust {rc:g}, rho_mantle {rm:g}")
    print('group | share of continent-throughout area | d_lith mean, p10/p50/p90 | |d_lith| > 100 / > 250 m share | '
          'airy mean | dyn mean | dz mean | gap mean | corr(d_lith, gap)')
    out = {}
    for name, m in groups.items():
        if m.sum() < 10:
            print(f'{name}: {m.sum()} parcels'); continue
        w = area[m]; dl = d_lith[m]
        r = np.corrcoef(dl, gap[m])[0, 1] if dl.std() > 0 else float('nan')
        row = {'share': float(w.sum()/area[cont_all].sum()), 'd_lith_mean': float(np.average(dl, weights=w)),
               'd_lith_q': wq(dl, w), 'gt100': float(w[np.abs(dl) > 100].sum()/w.sum()),
               'gt250': float(w[np.abs(dl) > 250].sum()/w.sum()), 'airy_mean': float(np.average(airy[m], weights=w)),
               'dyn_mean': float(np.average(dyn[m], weights=w)), 'dz_mean': float(np.average(dz[m], weights=w)),
               'gap_mean': float(np.average(gap[m], weights=w)), 'corr_gap': float(r),
               'abs_d_lith_mean': float(np.average(np.abs(dl), weights=w)),
               'abs_airy_mean': float(np.average(np.abs(airy[m]), weights=w)),
               'abs_dyn_mean': float(np.average(np.abs(dyn[m]), weights=w))}
        out[name] = row
        print(f"{name} | {row['share']:.3f} | {row['d_lith_mean']:+.0f}, " + ' / '.join(f'{v:+.0f}' for v in row['d_lith_q']) +
              f" | {row['gt100']:.3f} / {row['gt250']:.3f} | {row['airy_mean']:+.0f} | {row['dyn_mean']:+.0f} | "
              f"{row['dz_mean']:+.0f} | {row['gap_mean']:+.0f} | {r:+.2f}")
    print('mean |d_lith| against mean |airy| and mean |dyn| (the two channels the forcing carries):')
    for name, row in out.items():
        print(f"  {name}: {row['abs_d_lith_mean']:.0f} / {row['abs_airy_mean']:.0f} / {row['abs_dyn_mean']:.0f} m")
    print(json.dumps(out))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 30)
