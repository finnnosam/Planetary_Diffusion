"""Pass 55: why the planet's interior land (>= 600 km from the sea) stands ~740 m above Earth's. Terms, area-weighted
medians over interior land, six seeds, 0 Myr; the planet's own local Airy line (cont_ref 38 km at +300 m on the datum,
rho_c 2,830, rho_m 3,300) applied to crystalline crust, sediment at rho_s (bulk 2,400 as a round figure).
Earth: CRUST1.0 crsthk - sedthk (crystalline) and sedthk, ETOPO, ice masked; Earth's sea level is the datum."""
import sys, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import DATA, OUTPUTS
from planet_reviews.crust_review import run_params, sediment_bulk_m
from planet_reviews.coast_profile import distance_to_ocean_km
from planet_earth.reference import ice_mask
RC, RM, RS = 2830., 3300., 2400.
def wmed(v, w):
    o = np.argsort(v); c = np.cumsum(w[o]); return float(v[o][np.searchsorted(c, c[-1]/2)])
def airy(Tc, S): return 300. + (Tc-38000.)*(1-RC/RM) + S*(1-RS/RM)
g = CubeSphere(217); a = g.area
lat = np.degrees(np.arcsin(g.xyz[:, 2])); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ii, jj = np.clip((90-lat).astype(int), 0, 179), np.clip((lon+180).astype(int), 0, 359)
ld = lambda f: np.loadtxt(DATA/'crust1'/f)[:, 2].reshape(180, 360)[ii, jj]*1e3
T, S = ld('crsthk.xyz'), ld('sedthk.xyz')
ze = np.load(OUTPUTS/'planet-reference'/'earth-n217.npz')['elevation_m'].astype(float)
def report(name, z_sea, z_dat, Tc, S, land, far):
    m = far & land
    A = airy(Tc, S)
    r = dict(z_sea=wmed(z_sea[m], a[m]), z_datum=wmed(z_dat[m], a[m]), Tc=wmed(Tc[m], a[m])/1e3, S=wmed(S[m], a[m]),
             airy=wmed(A[m], a[m]), resid=wmed((z_dat-A)[m], a[m]), share=a[m].sum()/a[land].sum())
    print(f"{name:<16} z(sea) {r['z_sea']:6.0f}  z(datum) {r['z_datum']:6.0f}  crystalline {r['Tc']:5.1f} km  sediment {r['S']:5.0f} m"
          f"  Airy {r['airy']:6.0f}  z - Airy {r['resid']:6.0f}  (interior share of land {r['share']:.2f})")
    return r
keep = ~ice_mask(g); landE = (ze > 0) & keep
dE = distance_to_ocean_km(g, landE)
report('Earth', ze, ze, T-S, S, landE, dE >= 600)
for combo in sys.argv[2:]:
    rs = []
    for sd in sorted((Path(sys.argv[1])/combo).glob('seed*')):
        z = np.load(sd/'snapshots'/'t+0000.0.npz')
        e = z['elevation_m'].astype(float); sl = float(z['sea_level_m'])
        Sb = np.where(z['continental'], sediment_bulk_m(z), 0.); Tc = z['crust_thickness_m'].astype(float)
        land = e > 0; d = distance_to_ocean_km(g, land)
        rs.append(report(f'{combo} {sd.name}', e, e+sl, Tc, Sb, land, d >= 600))
    print(f"{combo} mean        " + '  '.join(f"{k} {np.mean([r[k] for r in rs]):.1f}" for k in rs[0]))
