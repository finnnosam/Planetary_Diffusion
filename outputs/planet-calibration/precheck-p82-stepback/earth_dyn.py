"""Pass 82 step-back, probe 8: does Earth's land sit on its dynamic-topography highs as the planet's does? Holdt et al.
2022's spherical-harmonic fits of observed residual topography (degree 1-30 and 1-10, 1-degree grids, km; oceanic
residual depths water-loaded, continental values from seismic crust: convention not confirmed, journal 73) sampled on
the n = 217 nodes: area means over the ocean, the land, and the land's elevation slices; the continental point
measurements beside them. The planet's field (dynamic_topography_m, six seeds) on the same slices. An indication only:
the continental part of Earth's field is the least constrained."""
import numpy as np, h5py, common as C
D = 'data/residual_topography/Residual depthtopographic measurements/Residual depthtopographic measurements/Holdt_et_al_2022_SupportingInformation/'
e = C.earth(); a = e['a']; z = e['z']; land = (z > 0) & ~e['ice']
def grid(name):
    with h5py.File(D+name, 'r') as f:
        keys = list(f.keys()); lon = f['lon'][:] if 'lon' in f else f['x'][:]; lat = f['lat'][:] if 'lat' in f else f['y'][:]
        v = f['z'][:]
    i = np.clip(np.round((e['lat']-lat[0])/(lat[1]-lat[0])).astype(int), 0, lat.size-1)
    j = np.round(((e['lon'] % 360)-lon[0])/(lon[1]-lon[0])).astype(int) % lon.size
    return v[i, j]*1e3, keys
SL = [('p20-30', .20, .30), ('p45-55', .45, .55), ('p85-95', .85, .95)]
for name in ('degree_1_to_30.grd', 'degree_1_to_10.grd'):
    f, keys = grid(name); ok = np.isfinite(f)
    print(f'{name} {keys}: global area mean {C.wmean(f[ok], a[ok]):.0f} m, rms {np.sqrt(C.wmean(f[ok]**2, a[ok])):.0f}; ocean (z < 0) mean {C.wmean(f[ok & (z < 0)], a[ok & (z < 0)]):.0f}; '
          f'land mean {C.wmean(f[ok & land], a[ok & land]):.0f}')
    q = {n: C.wq(z[land], a[land], [lo, hi]) for n, lo, hi in SL}
    print('   land slices (mean of the field, m): ' + ', '.join(f'{n} {C.wmean(f[land & ok & (z >= q[n][0]) & (z <= q[n][1])], a[land & ok & (z >= q[n][0]) & (z <= q[n][1])]):.0f}' for n, _, _ in SL))
pts = np.loadtxt(D+'continental_measurements_seismic_1_deg.xyz'); print(f'continental point measurements: {len(pts)} points, mean {1e3*pts[:, 2].mean():.0f} m, median {1e3*np.median(pts[:, 2]):.0f}, sd {1e3*pts[:, 2].std():.0f}')
o = []
for s in C.SEEDS:
    P = C.planet(s); zs = P['zd']-P['sea']; l = (zs > 0); ap = P['a']; d = P['dyn']; da = P['dyn_air']
    q = {n: C.wq(zs[l], ap[l], [lo, hi]) for n, lo, hi in SL}
    o.append([C.wmean(d, ap), np.sqrt(C.wmean(d*d, ap)), C.wmean(d[~l], ap[~l]), C.wmean(d[l], ap[l]), np.sqrt(C.wmean(da*da, ap))] + [C.wmean(d[l & (zs >= q[n][0]) & (zs <= q[n][1])], ap[l & (zs >= q[n][0]) & (zs <= q[n][1])]) for n, _, _ in SL])
o = np.array(o); m, s = o.mean(0), o.std(0, ddof=1)
print(f'planet (as applied): global mean {m[0]:.0f} m, rms {m[1]:.0f} (air-loaded rms {m[4]:.0f}); under the sea {m[2]:.0f}; land mean {m[3]:.0f} (sd seeds {s[3]:.0f}); land slices p20-30 {m[5]:.0f} ({s[5]:.0f}), p45-55 {m[6]:.0f} ({s[6]:.0f}), p85-95 {m[7]:.0f} ({s[7]:.0f})')
