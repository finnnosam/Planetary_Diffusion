"""Pass 69: continental lithosphere age at 0 Myr by distance from the coast, by seed; share reset during the run."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl

g = CubeSphere(217); a = g.area
BINS = [(0, 150), (150, 600), (600, 1200), (1200, 1e9)]
for arm in sorted(x for x in Path(sys.argv[1]).glob('*/') if (x/'seed1001').is_dir()):
    rows = []
    for sd in sorted(arm.glob('seed*')):
        z = np.load(sd/'snapshots'/'t+0000.0.npz')
        T = z['crust_thickness_m']+z['sediment_bulk_m']*z['continental']
        cont = z['continental'].astype(bool) & (T >= 20e3)
        land = z['elevation_m'] > 0
        d = gl.distance_to(~land, g.nbr, g.dist, 1e7)            # distance from the coast, on land
        age = z['lithosphere_age_myr'].astype(np.float64)
        out = []
        for lo, hi in BINS:
            m = cont & land & (d >= lo*1e3) & (d < hi*1e3) & np.isfinite(age)
            if not m.any():
                out.append('   --              '); continue
            o = np.argsort(age[m]); c = np.cumsum(a[m][o]); med = age[m][o][np.searchsorted(c, .5*c[-1])]
            young = a[m & (age < 250)].sum()/a[m].sum()
            out.append(f"{med:5.0f} ({young:.2f} <250, {a[m].sum()/a[cont].sum():.2f} of cont)")
        rows.append(f"  {sd.name}: " + " | ".join(out))
    print(arm.name, "  bins km: 0-150 | 150-600 | 600-1200 | 1200+ : median age Myr (share reset in the run, share of continent)")
    print("\n".join(rows), flush=True)
