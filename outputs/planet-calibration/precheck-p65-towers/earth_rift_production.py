"""Pass 65 (K): Earth's rifted-margin production since 250 Ma, on the n = 217 graph.

At each reconstruction time t (every 10 Myr), Mueller et al. (2019) paleo-age grid (0.1 deg, sampled at the nodes;
no age = not ocean = continent in the wide sense), the continental nodes beside ocean floor <= 10 Myr old, leaving out
those within 300 km of a subduction zone resolved from the same model at t (pygplates; their young floor is a
back-arc or a young plate arriving at a trench, not a breakup). Floor that young beside a margin at t was made in
(t, t+10] Ma, so length / 10 Myr is the production of rifted margin over that interval, each margin counted at its
birth (no survival bias). Both conjugate margins count (the instrument reads every continental node beside the new
ocean). Length = nodes x graph spacing, the same instrument as the planet reading (planet_reviews.rift_production).
Also the continent's area (non-ocean nodes) for the production per unit area.
  python earth_rift_production.py > earth_rift_production.txt   (writes earth_rift_production.json)"""
import json
import numpy as np
import h5py
import pygplates
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl

D = 'data/Muller2019/'
AG = D+'Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-%d.nc'
SF = D+'Muller_etal_2019_Tectonics_PlateMotionModel_v3.0/SimplifiedFiles/'
rot = pygplates.RotationModel(SF+'Muller_etal_2019_CombinedRotations.rot')
topo = pygplates.FeatureCollection(SF+'Muller_etal_2019_PlateBoundaries_DeformingNetworks.gpmlz')
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean()))
la = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lo = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
SZ = pygplates.FeatureType.gpml_subduction_zone


def nbv(m):
    return np.any(np.where(ok, m[np.where(ok, nb, 0)], False), 1)


def ages(t):
    f = h5py.File(AG % t)
    lon = f['lon'][:]; lat = f['lat'][:]; z = f['z']
    i = np.clip(np.round((la-lat[0])/(lat[1]-lat[0])).astype(int), 0, lat.size-1)
    j = np.clip(np.round((lo-lon[0])/(lon[1]-lon[0])).astype(int), 0, lon.size-1)
    a = z[:][i, j].astype(float)
    a[(a < 0) | (a > 1e3)] = np.nan
    return a


def subduction_nodes(t):
    rt, rs = [], []
    pygplates.resolve_topologies(topo, rot, rt, float(t), rs)
    pts = []
    for sec in rs:
        if sec.get_topological_section_feature().get_feature_type() != SZ:
            continue
        for sub in sec.get_shared_sub_segments():
            geom = sub.get_resolved_geometry()
            if geom is None:
                continue
            for p in geom.to_tessellated(np.radians(.2)).get_points():
                pts.append(p.to_xyz())
    m = np.zeros(g.size, bool)
    if pts:
        from planet_grid import fields
        m[fields.nearest_nodes(g, np.array(pts))] = True
    return m


rows = []
for t in range(0, 251, 10):
    a = ages(t)
    cont = ~np.isfinite(a)
    young = (~cont) & (a <= 10.)
    edge = cont & nbv(~cont)
    ye = cont & nbv(young)
    sz = subduction_nodes(t)
    near = gl.distance_to(sz, nb, g.dist, 400e3) <= 300e3 if sz.any() else np.zeros(g.size, bool)
    L_all = ye.sum()*sp/1e3; L_rift = (ye & ~near).sum()*sp/1e3
    A = g.area[cont].sum()/1e12
    rows.append({'t_ma': t, 'young_edge_km': L_all, 'rifted_km': L_rift, 'edge_km': edge.sum()*sp/1e3,
                 'cont_mkm2': A, 'prod_km_myr': L_rift/10., 'prod_per_mkm2': L_rift/10./A,
                 'subduction_nodes': int(sz.sum())})
    r = rows[-1]
    print(f"{t:4d} Ma: young edge {L_all/1e3:5.1f}k km, rifted (> 300 km from subduction) {L_rift/1e3:5.1f}k km -> "
          f"{r['prod_km_myr']:6.0f} km/Myr over ({t}, {t+10}] Ma; continent {A:5.1f} Mkm2 ({r['prod_per_mkm2']:.2f} km/Myr/Mkm2); "
          f"whole edge {r['edge_km']/1e3:5.1f}k; subduction nodes {r['subduction_nodes']}", flush=True)
P = np.array([r['prod_km_myr'] for r in rows[:25]]); Pa = np.array([r['prod_per_mkm2'] for r in rows[:25]])
print(f"mean 0-250 Ma: {P.mean():.0f} km/Myr ({Pa.mean():.3f} km/Myr/Mkm2); 250-150 Ma {P[15:].mean():.0f}, 150-50 {P[5:15].mean():.0f}, "
      f"50-0 {P[:5].mean():.0f}; made over 250 Myr {P.sum()*10/1e3:.0f}k km")
json.dump(rows, open('earth_rift_production.json', 'w'), indent=1)
