"""Pass 74: the planet's open-ocean sediment (bulk) by floor age at 0 Myr, same classes as earth_pelagic.py."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.ocean_sediment_probe import planet_fields, distances
g = CubeSphere(217)
for arm in sys.argv[1:]:
    rows = {}
    for sd in sorted(Path(arm).glob('seed*')):
        z, F = planet_fields(sd, 0, g)
        dt_, dc_ = distances(g, F)
        age = F['age']; sed = F['sed']
        for cmin, tmin in ((200, 200), (500, 300), (1200, 500)):
            m = F['ocean'] & (dc_ >= cmin) & (dt_ >= tmin)
            r = rows.setdefault((cmin, tmin), [])
            r.append((g.area[m].sum()/1e12, np.average(age[m], weights=g.area[m]), np.average(sed[m], weights=g.area[m]),
                      np.average(sed[F['ocean']], weights=g.area[F['ocean']]), (g.area*sed)[F['ocean']].sum()/1e15))
    for k, r in rows.items():
        r = np.array(r).mean(0)
        print(f'{Path(arm).name} continent >= {k[0]} km, trench >= {k[1]} km: area {r[0]:.1f} Mkm2, age {r[1]:.0f} Myr, sediment {r[2]:.0f} m bulk | whole ocean mean {r[3]:.0f} m, stock {r[4]:.0f} Mkm3 bulk')
