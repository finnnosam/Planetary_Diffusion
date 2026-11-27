# Pass 63: continental edge against young ocean floor through time: the length of continent beside ocean floor
# <= 10 Myr old, every 10 Myr (a breakup makes it; a ridge leaving the margin ends it). Earth: Mueller et al. (2019)
# age grids (0.1 degree, sampled at the n = 217 nodes; non-ocean = no age), planet: label-free continent (crust +
# sediment >= 20 km) beside ocean_age_myr <= 10. Also the whole continent-ocean edge for scale.
#   python young_edge.py [seed dirs...]
import sys, glob, numpy as np, h5py
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, run_params

AG = '/mnt/user-data/uploads/Planetary_Diffusion/data/Muller2019/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-%d.nc'
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean()))


def edges(cont, age):
    young = (~cont) & (age <= 10.)
    nbv = lambda m: np.any(np.where(ok, m[np.where(ok, nb, 0)], False), 1)
    return (cont & nbv(young)).sum()*sp/1e6, (cont & nbv(~cont)).sum()*sp/1e6


la = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lo = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
E = {}
for t in range(0, 251, 10):
    f = h5py.File(AG % t); lon = f['lon'][:]; lat = f['lat'][:]
    i = np.clip(np.round((la-lat[0])/(lat[1]-lat[0])).astype(int), 0, lat.size-1)
    j = np.clip(np.round((lo-lon[0])/(lon[1]-lon[0])).astype(int), 0, lon.size-1)
    age = f['z'][:][i, j].astype(float)
    E[-t] = edges(~np.isfinite(age), np.where(np.isfinite(age), age, np.inf))
P = []
for s in map(Path, sys.argv[1:]):
    p = run_params(s); r = {}
    for f in glob.glob(str(s/'snapshots'/'t*.npz')):
        z = np.load(f); t = float(z['time_myr'])
        if t % 10:
            continue
        cont = column_thickness_m(z, p) >= 20e3
        r[round(t)] = edges(cont, np.where(cont, np.inf, z['ocean_age_myr'].astype(float)))
    P.append((s, r))
print('t (Myr) | Earth: edge beside floor <= 10 Myr / whole edge (k km) | ' + ' | '.join(f'{s.parent.name}/{s.name}' for s, _ in P))
for t in range(-250, 1, 10):
    print('%5d   %6.1f / %6.1f   ' % (t, *E[t]) + '  '.join('%6.1f / %6.1f' % r[t] if t in r else '      -        ' for _, r in P))
ts = [t for t in range(-240, 1, 10)]
print('mean -240..0: Earth %.1f k km beside young floor (%.3f of its edge)' % (np.mean([E[t][0] for t in ts]), np.mean([E[t][0]/E[t][1] for t in ts])))
for s, r in P:
    tt = [t for t in ts if t in r]
    print('mean -240..0: %s %.1f k km (%.3f of its edge)' % (s, np.mean([r[t][0] for t in tt]), np.mean([r[t][0]/r[t][1] for t in tt])))
