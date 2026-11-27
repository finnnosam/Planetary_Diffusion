# Pass 67: where does the rifted margin go? For ALL continental edge (label-free, total >= 20 km, 46 km nodes), by
# margin age (oldest ocean floor within 300 km labelled to that edge node, as margin_age.py): length, share passive
# (> 300 km from any plate boundary), median distance to the boundary and to the ridge axis (ocean floor < 3 Myr,
# same rule on Earth and planet), and the implied half-spreading rate (ridge distance / margin age).
# Earth: CRUST1.0 + Bird PB2002 + Mueller et al. 2019 age grid at 0 Ma. Planet: snapshots.
#   python birth_probe.py [--t -40] <seed dirs...>     (no dirs: Earth only)
import sys, numpy as np, h5py
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import earth_nodes
from planet_earth.bird import earth_boundaries
from planet_tectonics.graphops import nearest_label
AGE = 'data/Muller2019/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc'
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean()))
BINS = [(0, 10), (10, 25), (25, 50), (50, 100), (100, 150), (150, 200), (200, 1e9)]

def probe(cont, bnd, age):
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    db = gl.distance_to(bnd, nb, g.dist, 5000e3)
    ridge = (~cont) & np.isfinite(age) & (age < 3)
    dr = gl.distance_to(ridge, nb, g.dist, 8000e3)
    lab, d = nearest_label(np.where(edge, np.arange(g.size), -1).astype(np.int64), nb, g.dist)
    sel = (~cont) & np.isfinite(age) & (d <= 300e3) & (lab >= 0)
    m = np.full(g.size, np.nan); np.fmax.at(m, lab[sel], age[sel])
    ma = m[edge]; dbe = db[edge]/1e3; dre = dr[edge]/1e3; pas = dbe > 300
    out = []
    for lo, hi in BINS + [(None, None)]:
        f = ~np.isfinite(ma) if lo is None else (ma >= lo) & (ma < hi)
        if not f.any(): out.append((lo, hi, 0, np.nan, np.nan, np.nan, np.nan)); continue
        rate = np.median(dre[f]/np.maximum(ma[f], 1)) if lo is not None else np.nan
        out.append((lo, hi, f.sum()*sp/1e6, pas[f].mean(), np.median(dbe[f]), np.median(dre[f]), rate/10))
    return out, edge.sum()*sp/1e6, pas.sum()*sp/1e6

def show(name, res):
    out, tot, pas = res
    print(f'{name}: edge {tot:.1f}k km, passive {pas:.1f}k ({pas/tot:.2f})')
    for lo, hi, L, fp, mb, mr, v in out:
        lab = 'undated' if lo is None else f'{lo:.0f}-{hi:.0f}' if hi < 1e8 else f'{lo:.0f}+'
        print(f'   {lab:>8s} Myr: {L:5.1f}k km, passive {fp:4.2f}, to boundary p50 {mb:5.0f} km, to ridge p50 {mr:5.0f} km, half rate {v:4.1f} cm/yr')

def mean(rs):
    o = [rs[0][0][i][:2] + tuple(np.nanmean([r[0][i][j] for r in rs]) for j in range(2, 7)) for i in range(len(BINS)+1)]
    return o, np.mean([r[1] for r in rs]), np.mean([r[2] for r in rs])

if __name__ == '__main__':
    a = sys.argv[1:]; t = 0.
    if a[:1] == ['--t']: t = float(a[1]); a = a[2:]
    if not a:
        T, _ = earth_nodes(g); cont = label_free_continent(T)
        f = h5py.File(AGE); lon = f['lon'][:]; lat = f['lat'][:]; z = f['z'][:]
        la = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lo = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
        i = np.clip(np.round((la-lat[0])/(lat[1]-lat[0])).astype(int), 0, lat.size-1)
        j = np.clip(np.round((lo-lon[0])/(lon[1]-lon[0])).astype(int), 0, lon.size-1)
        age = z[i, j].astype(float); age[cont] = np.nan
        show('Earth (CRUST1.0, Bird, Mueller 2019)', probe(cont, earth_boundaries(g), age))
    R = []
    for s in map(Path, a):
        z = np.load(s/'snapshots'/f't{t:+07.1f}.npz'); p = run_params(s)
        Tcol = column_thickness_m(z, p); bnd = z['boundary_type'] != 0
        cont = label_free_continent(Tcol)
        age = np.where(cont, np.nan, z['ocean_age_myr'].astype(float))
        R.append(probe(cont, bnd, age))
    if R: show(f'{Path(a[0]).parent.name} ({len(R)} seeds, t {t:+.0f} Myr, means)', mean(R))
