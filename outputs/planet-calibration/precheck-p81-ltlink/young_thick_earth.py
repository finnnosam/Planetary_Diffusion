"""Pass 81: Earth's seismic L (SL2013sv) under Cenozoic-class crust (CRUST1.0 types of age class 33 Ma: orogens,
arcs, rifts, extended crust) by total crust bin and region. Run from the repo root."""
import numpy as np
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_earth.reference import ice_mask
from planet_tectonics import lithosphere as lt
g = CubeSphere(217)
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]
tn = np.array([l.split() for l in (DATA/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
cont = ~np.isin(tn, lt.OCEANIC_TYPES) & (T >= 20) & ~ice_mask(g)
dd = np.loadtxt(DATA.joinpath(*lt.HOGGARD_SL2013SV), delimiter=',')
grid = np.full((361, 721), np.nan)
grid[np.round((90-dd[:, 1])*2).astype(int), np.round((dd[:, 0] % 360)*2).astype(int) % 720] = dd[:, 2]
grid[:, 720] = grid[:, 0]
L = grid[np.round((90-lat)*2).astype(int), np.round((lon % 360)*2).astype(int)]
age = np.array([lt.CLASS_AGE_MA.get(t, np.nan) for t in tn])
young = cont & (age == 33.) & np.isfinite(L)
regions = {'Andes (lon -82..-60)': (lon > -82) & (lon < -60), 'Tibet-Himalaya (lon 70..105, lat 25..45)':
           (lon > 70) & (lon < 105) & (lat > 25) & (lat < 45), 'Iran-Anatolia-Caucasus (lon 25..70)': (lon > 25) & (lon < 70),
           'N America Cordillera (lon -130..-100)': (lon > -130) & (lon < -100), 'all': np.ones(g.size, bool)}
types = sorted(set(tn[young]))
print('Cenozoic-class types present:', types)
for name, rg in regions.items():
    out = []
    for lo, hi in [(20, 35), (35, 45), (45, 55), (55, 80)]:
        k = young & rg & (T >= lo) & (T < hi)
        out.append(f'{lo}-{hi}: {np.average(L[k], weights=g.area[k]) if k.any() else float("nan"):.0f} km ({g.area[k].sum()/1e12:.2f} Mkm2)')
    print(f'{name}: ' + '; '.join(out))
for t in types:
    k = young & (tn == t)
    kt = k & (T >= 45)
    print(f'  type {t}: area {g.area[k].sum()/1e12:.2f} Mkm2, L {np.average(L[k], weights=g.area[k]):.0f}; crust >= 45: {g.area[kt].sum()/1e12:.2f} Mkm2, L {np.average(L[kt], weights=g.area[kt]) if kt.any() else float("nan"):.0f}')
