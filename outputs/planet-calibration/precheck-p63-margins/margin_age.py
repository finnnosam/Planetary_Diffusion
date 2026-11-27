# Pass 63: how old are the passive margins? Earth (CRUST1.0 label-free edge on n = 217, Bird PB2002 boundaries,
# Mueller et al. 2019 age grid at 0 Ma) against planet runs at 0 Myr (label-free, 46 km and the 1 degree instrument).
# Margin age of a passive-edge node = the oldest ocean floor among the oceanic nodes within 300 km whose nearest
# passive-edge node it is (the floor next to a rifted margin dates its breakup). Also: rift events in the run.
#   python margin_age.py <seed dirs...>
import sys, json, numpy as np, h5py
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import earth_nodes, to_crust1_cells, BOUNDARY_KM
from planet_earth.bird import earth_boundaries
from planet_tectonics.graphops import nearest_label

AGE = '/mnt/user-data/uploads/Planetary_Diffusion/data/Muller2019/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-0.nc'
g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))


def margin_age(cont, bnd, age):
    edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
    far = gl.distance_to(bnd, nb, g.dist, 5000e3) > BOUNDARY_KM*1e3
    pas = edge & far
    lab0 = np.where(pas, np.arange(g.size), -1).astype(np.int64)
    lab, d = nearest_label(lab0, nb, g.dist)
    sel = (~cont) & np.isfinite(age) & (d <= 300e3) & (lab >= 0)
    m = np.full(g.size, np.nan)
    np.fmax.at(m, lab[sel], age[sel])
    ma = m[pas]
    f = np.isfinite(ma)
    q = np.percentile(ma[f], [25, 50, 75]) if f.any() else [np.nan]*3
    return dict(passive_k_km=pas.sum()*sp/1e6, active_k_km=(edge & ~far).sum()*sp/1e6, cont_mkm2=a[cont].sum()/1e12,
                dated=f.mean(), mean=np.nanmean(ma), p25=q[0], p50=q[1], p75=q[2],
                young50=(ma[f] < 50).mean(), old100=(ma[f] >= 100).mean(),
                passive_per_mkm2=pas.sum()*sp/1e3/(a[cont].sum()/1e12),
                **{f'L{lo}': pas.sum()*sp/1e6*((ma[f] >= lo) & (ma[f] < lo+50)).mean() for lo in (0, 50, 100, 150, 200)})


def fmt(name, r):
    return (f"{name:28s} passive {r['passive_k_km']:6.1f}k km  active {r['active_k_km']:6.1f}k  cont {r['cont_mkm2']:6.1f} Mkm2  "
            f"passive per Mkm2 {r['passive_per_mkm2']:5.0f} km  margin age mean {r['mean']:5.1f} p25/50/75 {r['p25']:4.0f} / "
            f"{r['p50']:4.0f} / {r['p75']:4.0f} Myr  <50 Myr {r['young50']:.2f}  >=100 Myr {r['old100']:.2f}  dated {r['dated']:.2f}"
            + '\n      passive length by margin age 0-50/50-100/100-150/150-200/200-250 Myr (k km): ' + ' / '.join('%.1f' % r[f'L{lo}'] for lo in (0, 50, 100, 150, 200)))


def earth():
    T, _ = earth_nodes(g)
    cont = label_free_continent(T)
    f = h5py.File(AGE); lon = f['lon'][:]; lat = f['lat'][:]; z = f['z'][:]
    la = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lo = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    i = np.clip(np.round((la-lat[0])/(lat[1]-lat[0])).astype(int), 0, lat.size-1)
    j = np.clip(np.round((lo-lon[0])/(lon[1]-lon[0])).astype(int), 0, lon.size-1)
    age = z[i, j].astype(float)
    return margin_age(cont, earth_boundaries(g), age)


if __name__ == '__main__':
    print(fmt('Earth (CRUST1.0, 46 km)', earth()))
    R46, R1 = [], []
    for s in map(Path, sys.argv[1:]):
        z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s)
        Tcol = column_thickness_m(z, p); bnd = z['boundary_type'] != 0
        age = np.where(z['continental'].astype(bool), np.nan, z['ocean_age_myr'].astype(float))
        r46 = margin_age(label_free_continent(Tcol), bnd, age)
        T1 = to_crust1_cells(g, Tcol, z['continental'].astype(bool))[0]
        r1 = margin_age(label_free_continent(T1), bnd, age)
        ev = json.load(open(s/'stats.json')).get('event_list', [])
        kinds = {}
        for e in ev: kinds[e['kind']] = kinds.get(e['kind'], 0)+1
        print(fmt(f'{s.parent.name}/{s.name} 46 km', r46)); print(fmt('   1 degree', r1)); print('   events', kinds)
        R46.append(r46); R1.append(r1)
    if R46:
        mean = lambda R: {k: float(np.mean([r[k] for r in R])) for k in R[0]}
        print(fmt('mean, 46 km', mean(R46))); print(fmt('mean, 1 degree', mean(R1)))
