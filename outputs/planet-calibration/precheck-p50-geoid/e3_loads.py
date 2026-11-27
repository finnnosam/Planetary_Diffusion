"""Pass 50: E3's record (Müller 2019 through D's unchanged SlabRecord; same steps as precheck-p49-e3/e3_run.py) with the
present-day load layers saved (surface-density alm per 50 km layer, no alpha factor), so any viscosity profile's
kernels can be applied to the same loads. Also saves the loads at 10, 20 and 50 Ma."""
import sys
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_tectonics.params import PARAMS as P

R = 6371e3
T = np.load(sys.argv[1]); A = np.load(sys.argv[2]); out = sys.argv[3]
age = A['age_myr']; assert age.shape == T['t'].shape
g = CubeSphere(217)
z_cut = P['flow_z_cut_km']*1e3; lmax = int(P['flow_lmax'])
radii = mf.load_radii(R, z_cut)
rec = mf.SlabRecord(P['slab_dip_deg'], P['lm_sink_mm_yr'])
tt = T['t']; save = {}
for t in range(250, -1, -1):
    rec.advance(1., R)
    m = (tt == t) & np.isfinite(age) & (T['conv_mm_yr'] > 0)
    conv = T['conv_mm_yr'][m]*1e-3
    mass = mf.thermal_mass_excess(age[m], P['rho_mantle'], P['rho_water'])*conv*1e6*T['length_km'][m]*1e3
    rec.deposit(T['xyz'][m], T['n_over'][m], mass, conv, -float(t))
    rec.merge(g, z_cut)
    if t in (50, 20, 10, 0):
        save[f'loads_{t}'] = rec.layer_loads(z_cut, len(radii), R, lmax, 1, None)
        save[f'mass_by_layer_{t}'] = rec.mass_by_layer(z_cut, len(radii))
        print(t, rec.total_mass(), flush=True)
np.savez_compressed(out, radii=radii, lmax=lmax, deposited=rec.deposited_kg, exited=rec.exited_kg, held=rec.total_mass(), **save)
print('done')
