"""Pass 111 probe: oceanic crust thinner than ocean_thick_m (7 km) in saved snapshots, which no ledger books (the
invariant of tests/test_plates.py::test_short_run_invariants and test_shelf_process_keeps_the_volume_ledger).
Usage (repo root): python ocean_thin.py <arm folder> [<arm folder> ...]
Prints per arm and snapshot (-200 and 0 Myr, six seeds): the minimum oceanic thickness, the area share of the ocean
below 6,999 m and below 6,900 m, and the deficit below 7 km in km^3 (seed mean); then for seed 1001 every listed
snapshot: thinned nodes, deficit, share next to continental crust, boundary types (0 interior), median ocean age."""
import sys
import numpy as np
from planet_grid.sphere_graph import CubeSphere

g = CubeSphere(217); a = g.area
for arm in sys.argv[1:]:
    for t in ('t-0200.0', 't+0000.0'):
        r = []
        for s in range(1001, 1007):
            d = np.load(f'{arm}/seed{s}/snapshots/{t}.npz'); oc = ~d['continental']
            T = d['crust_thickness_m'][oc].astype(float); ao = a[oc]
            r.append((T.min(), ao[T < 6999].sum()/ao.sum(), ao[T < 6900].sum()/ao.sum(),
                      (np.clip(7000-T, 0, None)*ao).sum()/1e9))
        r = np.array(r)
        print(f'{arm} {t}: min {r[:, 0].min():.0f} m | share < 6999 {r[:, 1].mean():.4f} | < 6900 {r[:, 2].mean():.5f}'
              f' | deficit {r[:, 3].mean():.0f} km3')
    for t in ('t-0240.0', 't-0230.0', 't-0220.0', 't-0100.0', 't+0000.0'):
        d = np.load(f'{arm}/seed1001/snapshots/{t}.npz')
        c = d['continental']; T = d['crust_thickness_m']; thin = ~c & (T < 6999)
        if not thin.any():
            print(f'  seed1001 {t}: none'); continue
        nb = g.nbr[thin]; nearc = np.where(nb >= 0, c[np.maximum(nb, 0)], False).any(1)
        print(f'  seed1001 {t}: {thin.sum()} nodes, deficit {(np.clip(7000-T[thin], 0, None)*a[thin]).sum()/1e9:.0f} km3,'
              f' next to continent {nearc.mean():.2f}, boundary types {np.bincount(d["boundary_type"][thin], minlength=5)},'
              f' ocean age p50 {np.nanmedian(d["ocean_age_myr"][thin]):.0f} Myr')
