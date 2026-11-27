# Pass 60: Earth's thin crust of oceanic types (A0 A1 V1 Y3, B-, W-) at total >= 20 km: where, how much is sediment
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_earth.crust1 import C1
from planet_earth.reference import ice_mask
g = CubeSphere(217); a = g.area/1e12
e = np.load('outputs/planet-calibration/precheck-p60-area/earth_thin_nodes.npz'); T, ty, z = e['T'], e['ty'], e['z']
S = np.loadtxt(C1/'sedthk.xyz')[:, 2].reshape(180, 360)*1e3
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359); S = S[r, c]
ice = ice_mask(g)
thin = (T >= 20e3) & (T < 30e3)
REG = [('Antarctica/Greenland ice', ice), ('Arctic >66N', lat > 66), ('Bay of Bengal', (lat > 0) & (lat < 23) & (lon > 80) & (lon < 95)),
       ('Gulf of Mexico/Caribbean', (lat > 8) & (lat < 31) & (lon > -98) & (lon < -60)), ('W Pacific 140E-180', (lon > 140) & (lat > -40) & (lat < 45)),
       ('Indian Ocean', (lon > 40) & (lon < 110) & (lat < 0) & (lat > -60)), ('S Atlantic', (lon > -60) & (lon < 20) & (lat < 0) & (lat > -60))]
for grp, codes in (('oceanic A0 A1 V1 Y3', ['A0', 'A1', 'V1', 'Y3']), ('plateau B- W-', ['B-', 'W-']), ('shelf C-', ['C-']), ('slope S- T- U-', ['S-', 'T-', 'U-'])):
    m = thin & np.isin(ty, codes)
    w = a[m]; print('%s: %.1f Mkm2, sediment p50 %.1f km, crystalline p50 %.1f km, lat p50 %.0f' % (grp, w.sum(), np.median(S[m])/1e3, np.median((T-S)[m])/1e3, np.median(lat[m])))
    left = m.copy()
    for name, rm in REG:
        k = left & rm; left &= ~rm
        if a[k].sum() > .2: print('    %-26s %.1f' % (name, a[k].sum()))
    print('    %-26s %.1f' % ('elsewhere', a[left].sum()))
m = thin & ~ice
print('thin crust outside the ice regions: %.1f Mkm2 (continent outside ice %.1f; share %.3f)' % (a[m].sum(), a[(T >= 20e3) & ~ice].sum(), a[m].sum()/a[(T >= 20e3) & ~ice].sum()))
