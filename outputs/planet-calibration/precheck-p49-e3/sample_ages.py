"""E3 step 2 (pass 49), runs where the grids are: the subducting floor's age at each trench point from Müller et al.
(2019) v2.0 age grids (Muller_etal_2019_Tectonics_v2.0_AgeGrid-<t>.nc, netCDF-4 read with h5py, or netCDF-3 by a
minimal reader). points.npz: t (N), lon (3, N), lat (3, N) at 0.5 / 1.0 / 1.5 deg from the trench into the
subducting plate. Nearest grid node; the first finite of the three offsets. Writes ages_part_<t0>_<t1>.npz.

  python3 sample_ages.py points.npz <grid folder> <t0> <t1> <out folder>
"""
import sys
from pathlib import Path
import numpy as np


def read_grid(path):
    import h5py
    with h5py.File(path, 'r') as f:
        names = list(f.keys())
        xn = next(n for n in names if n.lower() in ('lon', 'longitude', 'x'))
        yn = next(n for n in names if n.lower() in ('lat', 'latitude', 'y'))
        zn = next(n for n in names if n not in (xn, yn) and f[n].ndim == 2)
        z = f[zn][...].astype(np.float64)
        fill = f[zn].attrs.get('_FillValue')
        if fill is not None:
            z[z == np.asarray(fill).ravel()[0]] = np.nan
        return f[xn][...].astype(np.float64), f[yn][...].astype(np.float64), z


def sample(lon, lat, x, y, z):
    lon = np.where(lon < x.min()-1e-9, lon+360, lon); lon = np.where(lon > x.max()+1e-9, lon-360, lon)
    i = np.clip(np.rint((lat-y[0])/(y[1]-y[0])).astype(int), 0, len(y)-1)
    j = np.clip(np.rint((lon-x[0])/(x[1]-x[0])).astype(int), 0, len(x)-1)
    return z[i, j]


def main(points, folder, t0, t1, outdir):
    P = np.load(points)
    t, LON, LAT = P['t'], P['lon'], P['lat']
    folder = Path(folder)
    sel = (t >= t0) & (t <= t1)
    idx = np.flatnonzero(sel)
    age = np.full(idx.size, np.nan); which = np.full(idx.size, -1)
    for tt in range(int(t0), int(t1)+1):
        f = next(folder.rglob(f'*AgeGrid-{tt}.nc'))
        x, y, z = read_grid(f)
        m = t[idx] == tt
        k = idx[m]
        a = np.full(k.size, np.nan); w = np.full(k.size, -1)
        for o in range(LON.shape[0]):
            v = sample(LON[o, k], LAT[o, k], x, y, z)
            new = np.isnan(a) & np.isfinite(v)
            a[new] = v[new]; w[new] = o
        age[m] = a; which[m] = w
    np.savez_compressed(Path(outdir)/f'ages_part_{int(t0):03d}_{int(t1):03d}.npz', idx=idx, age_myr=age, offset=which)
    print(f'{t0}-{t1}: {idx.size} points, {np.isfinite(age).mean():.3f} with an age')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
