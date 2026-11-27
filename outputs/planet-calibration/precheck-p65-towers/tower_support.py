"""Pass 65: how much of each run's highest node is local Airy (crust + sediment) and how much is the elastic plate
holding up a narrow load (z - Airy - dyn). Airy as surface.elevation's own_airy (cont_ref 38 km at +300 m,
rho 2830 / 3300), continental sediment added at its bulk height x (1 - 2400/3300) as an Airy load; sea_level_m is
the sea on the datum, elevation_m is relative to it. Also over all nodes above 6 km: mean support."""
import sys
from pathlib import Path
import numpy as np
rc, rm, rs = 2830., 3300., 2400.
sw = Path(sys.argv[1])
for c in sorted(sw.glob('combo-0*')):
    rows=[]
    for s in sorted(c.glob('seed*')):
        best=None
        for f in sorted((s/'snapshots').glob('t*.npz')):
            z=np.load(f); m=float(z['elevation_m'].max())
            if best is None or m>best[0]: best=(m,f)
        z=np.load(best[1]); e=z['elevation_m'].astype(float)
        T=z['crust_thickness_m'].astype(float); S=z['sediment_bulk_m'].astype(float); d=z['dynamic_topography_m'].astype(float)
        airy=300.+(T-38000.)*(rm-rc)/rm+S*(1-rs/rm)-float(z['sea_level_m'])
        sup=e-airy-d
        i=int(e.argmax()); hi=(e>6000)&z['continental']
        rows.append((s.name, e[i], airy[i], d[i], sup[i], sup[hi].mean() if hi.any() else np.nan, np.percentile(sup[z['continental']&(e>0)],99)))
    print(f"== {c.name}")
    for r in rows: print("  %s top %.0f = Airy %.0f + dyn %.0f + plate %.0f | nodes > 6 km: plate mean %.0f | land plate p99 %.0f" % r)
