"""Pass 55: how high does continental crust stand at a given total thickness, planet vs Earth (n = 217, 0 Myr)?
Elevation relative to sea level, label-free continent (total >= 20 km), weighted medians per thickness bin.
Planet columns split into: sea on the datum, dynamic topography, sediment share. Earth: ETOPO (ice masked) on CRUST1.0."""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import DATA, OUTPUTS
from planet_reviews.crust_review import run_params, total_thickness_m, sediment_bulk_m
from planet_earth.reference import ice_mask

def wmed(v, w):
    if len(v) == 0: return np.nan
    o = np.argsort(v); c = np.cumsum(w[o]); return float(v[o][np.searchsorted(c, c[-1]/2)])

BINS = [20, 25, 30, 33, 36, 39, 42, 45, 50, 60]
g = CubeSphere(217); a = g.area
lat = np.degrees(np.arcsin(g.xyz[:, 2])); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ze = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360)
S = sum(np.loadtxt(DATA/'crust1'/f'sedthk.xyz')[:, 2].reshape(180, 360) for _ in [0]) if (DATA/'crust1'/'sedthk.xyz').exists() else None
ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
Te = T[ii, jj]; keep = ~ice_mask(g)
print('total km    ', ' '.join(f'{lo:>2}-{hi:<3}' for lo, hi in zip(BINS[:-1], BINS[1:])))
def row(name, z, Tk, m, extra=None):
    out = []
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        s = m & (Tk >= lo) & (Tk < hi); out.append(wmed(z[s], a[s]))
    print(f'{name:<12}', ' '.join(f'{v:6.0f}' for v in out))
    return out
earth = row('Earth z', ze, Te, keep & (Te >= 20))
if S is not None: row('Earth sed', S[ii, jj]*1e3, Te, keep & (Te >= 20))
sw = Path(sys.argv[1])
for combo in sys.argv[2:]:
    zs, ss, ds, sl, sd, cr = [], [], [], [], [], []
    for sd_ in sorted((sw/combo).glob('seed*')):
        z = np.load(sd_/'snapshots'/'t+0000.0.npz'); p = run_params(sd_)
        Tt = total_thickness_m(z, p)/1e3; m = Tt >= 20
        zs.append(row(f'{combo[-2:]} {sd_.name[-4:]} z', z['elevation_m'].astype(float), Tt, m))
    print(f'{combo} mean z ', ' '.join(f'{v:6.0f}' for v in np.nanmean(zs, 0)), '| minus Earth', ' '.join(f'{v:5.0f}' for v in np.nanmean(zs, 0)-earth))
    for sd_ in sorted((sw/combo).glob('seed*')):
        z = np.load(sd_/'snapshots'/'t+0000.0.npz'); p = run_params(sd_)
        Tt = total_thickness_m(z, p)/1e3; m = Tt >= 20
        sl.append(float(z['sea_level_m'])); sd.append(row('   dyn', z['dynamic_topography_m'].astype(float), Tt, m))
        ss.append(row('   sed bulk', np.where(z['continental'], sediment_bulk_m(z), 0.), Tt, m))
    print(f'{combo} sea on datum mean {np.mean(sl):.0f} m; dyn by bin', ' '.join(f'{v:5.0f}' for v in np.nanmean(sd, 0)), '; sed bulk', ' '.join(f'{v:5.0f}' for v in np.nanmean(ss, 0)))
