"""Pass 83 (journal 108): where the sediment lies, Earth against the planet, on the n = 217 nodes.
Earth: CRUST1.0 total (crust + sediment) and sediment (bulk) on the continent; GlobSed v3 (bulk) on the floor where
valid, else CRUST1.0; ETOPO for above / below sea; Seton 2020 floor ages. Label-free continent: total >= 20 km.
Planet: a sweep arm's 0 Myr snapshots (sediment_bulk_m; continent = crust + bulk sediment >= 20 km).
Domains: continent above sea; continent under the sea; off the continent by distance from it (0-200, 200-500,
500-1000, > 1000 km); and the far floor by age. Volumes in Mkm^3 bulk, thickness m (area means).
Usage (repo root): python sediment_where.py <arm> [<arm> ...]"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl

g = CubeSphere(217); a = g.area; R = 6371e3
AREA = 4*np.pi*R**2
DBINS = [(0, 200), (200, 500), (500, 1000), (1000, 1e9)]
ABINS = [(0, 10), (10, 30), (30, 60), (60, 100), (100, 400)]

def earth():
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
    c1 = DATA/'crust1'
    T = np.loadtxt(c1/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
    S = np.loadtxt(c1/'sedthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
    Rf = OUTPUTS/'planet-reference'
    z = np.load(Rf/'earth-n217.npz')['elevation_m'].astype(float)
    ag = np.load(Rf/'earth-seafloor-age-n217.npz'); age = ag['age_myr'].astype(float)
    gs = np.load(Rf/'globsed-n217.npz'); h = gs['sediment_m'].astype(float)
    cont = T >= 20e3
    sed = np.where(~cont & np.isfinite(h) & (gs['valid_share'] >= .5), h, S)
    age = np.where(ag['dated_share'] >= .99, age, np.nan)
    return dict(z=z, cont=cont, sed=sed, age=age, sea=0.)

def planet(arm, seed):
    d = np.load(f'{arm}/seed{seed}/snapshots/t+0000.0.npz')
    sed = d['sediment_bulk_m'].astype(float); T = d['crust_thickness_m'].astype(float)
    cont = (T+sed) >= 20e3
    oa = d['ocean_age_myr'].astype(float)
    return dict(z=d['elevation_m'].astype(float), cont=cont, sed=sed, age=np.where(cont, np.nan, oa), sea=0.)

def table(e):
    dist = gl.distance_to(e['cont'], g.nbr, g.dist, 3e6)/1e3
    out = {}
    def put(k, m):
        A = a[m].sum(); V = (e['sed'][m]*a[m]).sum()
        out[k] = (A/a.sum(), V*AREA/a.sum()/1e9/1e6, V/A if A else np.nan)   # share, Mkm^3, mean m
    put('continent above sea', e['cont'] & (e['z'] >= 0))
    put('continent under sea', e['cont'] & (e['z'] < 0))
    for lo, hi in DBINS:
        put(f'floor {lo:.0f}-{hi:.0f} km' if hi < 1e8 else f'floor > {lo:.0f} km', ~e['cont'] & (dist >= lo) & (dist < hi))
    far = ~e['cont'] & (dist >= 500)
    for lo, hi in ABINS:
        put(f'  far floor {lo}-{hi} Myr', far & (e['age'] >= lo) & (e['age'] < hi))
    near = ~e['cont'] & (dist < 500)
    for lo, hi in ABINS:
        put(f'  near floor {lo}-{hi} Myr', near & (e['age'] >= lo) & (e['age'] < hi))
    put('all', np.ones(a.size, bool))
    return out

if __name__ == '__main__':
    E = table(earth())
    cols = [('Earth', E)]
    for arm in sys.argv[1:]:
        ts = [table(planet(arm, s)) for s in range(1001, 1007)]
        m = {k: tuple(np.nanmean([t[k][i] for t in ts]) for i in range(3)) for k in ts[0]}
        cols.append((Path(arm).name, m))
    print('domain | ' + ' | '.join(f'{n}: share / Mkm3 / m' for n, _ in cols))
    for k in E:
        print(f'{k:24s} | ' + ' | '.join(f'{c[k][0]:.3f} / {c[k][1]:6.1f} / {c[k][2]:6.0f}' for _, c in cols))
