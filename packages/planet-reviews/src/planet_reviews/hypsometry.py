"""Hypsometry and shelf measures shared by the sweeps (scripts/sweep_plates.py) and the elevation check."""
import numpy as np

QS = ('25', '50', '75', '90')
# rock delivered from the continents to the ocean, km^3/yr: ~14 Gt/yr of pre-human suspended load plus ~4 Gt/yr
# dissolved (Milliman & Farnsworth 2011; Syvitski et al. 2005), at 2.65 t/m^3, ~6.8; band kept wide
EARTH_EXPORT_KM3_YR = (4., 8.)


def shelf_metrics(z, cont, g):
    a = g.area
    land = z > 0
    return {'whole_p70_m': float(np.percentile(z, 70)), 'continental_submerged': float(a[cont & ~land].sum()/a[cont].sum()),
            'land_p25_m': float(np.percentile(z[land], 25)), 'shallow_sea_share': float(a[(z > -200) & ~land].sum()/a.sum()),
            'continental_fraction': float(a[cont].sum()/a.sum())}
