"""Pass 65: what the towers of sweep-20261104-214903 are made of. For each combo and seed, at the snapshot holding
the run's highest node: the top node's elevation, crust, sediment (bulk), dyn topo, boundary type, orogeny age,
and the area above 8,850 m / 6,000 m, plus the crust and sediment p99 on continental nodes."""
import sys, glob
from pathlib import Path
import numpy as np
sw = Path(sys.argv[1])
for c in sorted(sw.glob('combo-0*')):
    print(f"== {c.name}")
    for s in sorted(c.glob('seed*')):
        best=None
        for f in sorted((s/'snapshots').glob('t*.npz')):
            z=np.load(f); m=float(z['elevation_m'].max())
            if best is None or m>best[0]: best=(m,f)
        z=np.load(best[1]); e=z['elevation_m']; i=int(e.argmax())
        n=e.size; a=510.1e6/n  # km2 per node (equal-area assumed)
        cont=z['continental']
        T=z['crust_thickness_m']; S=z['sediment_bulk_m']
        hi=e>8850; mid=e>6000
        print(f"  {s.name} t={float(z['time_myr']):+.0f}: top {e[i]:.0f} m crust {T[i]/1e3:.1f} km sed {S[i]/1e3:.1f} km dyn {z['dynamic_topography_m'][i]:.0f} m bt {z['boundary_type'][i]} orog {z['orogeny_age_myr'][i]:.0f} thick+ {z['tectonic_thickening_m'][i]:.0f} m | >8.85 km {hi.sum()*a/1e3:.0f}k km2 (crust p50 {np.median(T[hi])/1e3 if hi.any() else 0:.0f} sed p50 {np.median(S[hi])/1e3 if hi.any() else 0:.1f}) >6 km {mid.sum()*a/1e3:.0f}k km2 | cont crust p99 {np.percentile(T[cont],99)/1e3:.1f} max {T[cont].max()/1e3:.1f} sed p99 {np.percentile(S[cont],99)/1e3:.1f} max {S.max()/1e3:.1f} km")
