"""z - Airy by band, split: C term, dynamic topography, sediment lighter than 2,400 (bulk), the rest (flexure, water)."""
import sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import coast_budget as cb
from planet_reviews.crust_review import run_params, sediment_bulk_m
from planet_reviews.coast_profile import BINS_KM, distance_to_ocean_km
from planet_tectonics.lithosphere import elevation_term
g = CubeSphere(217); A = g.area
for sd in sys.argv[1:]:
    p = run_params(sd)
    z = np.load(f'{sd}/snapshots/t+0000.0.npz'); c = z['continental'].astype(bool)
    e = z['elevation_m'].astype(float); zd = e+float(z['sea_level_m']); land = e > 0
    S = np.where(c, sediment_bulk_m(z), 0.); Ss = np.where(c, z['sediment_m'].astype(float), 0.)
    mass = p['rho_grain']*Ss+p['rho_water']*(S-Ss)
    sedlight = (p['rho_sediment']*S-mass)/p['rho_mantle']
    Tc = z['crust_thickness_m'].astype(float)
    res = zd-cb.airy_m(Tc, S, p)
    cterm = np.where(c, elevation_term(z['lithosphere_age_myr'].astype(float), p), 0.) if p.get('lithosphere') == 'thermal_age' else 0*res
    dyn = z['dynamic_topography_m'].astype(float)
    rest = res-cterm-dyn-sedlight
    d = distance_to_ocean_km(g, land)
    print(sd)
    for lo, hi in zip(BINS_KM[:-1], BINS_KM[1:]):
        m = land & (d >= lo) & (d < hi); w = A[m]
        f = lambda v: cb.wmed(v[m], w)
        mn = lambda v: float((v[m]*w).sum()/w.sum())
        print(f'  {lo:5.0f}-{hi:<6.0f} z-Airy p50 {f(res):6.0f} mean {mn(res):6.0f} | means: C {mn(cterm):6.0f} dyn {mn(dyn):6.0f} light sediment {mn(sedlight):6.0f} rest {mn(rest):6.0f}')
