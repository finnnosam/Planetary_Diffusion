"""E3 step 3 (pass 49): Earth's subduction history (Müller et al. 2019 v3.0) through D's own slab record and kernels.

Uses planet_tectonics.mantle_flow unchanged: SlabRecord (dip 58 deg toward the overriding side at the convergence rate,
12 mm/yr below 660 km, out at 2,850 km), merge on CubeSphere(217) nodes, thermal_mass_excess (the floor's GDH1
subsidence x (rho_m - rho_w)), alpha_ratio, FlowKernels with data/mantle/visc.sh08, cut 300 km, degrees 1-31 —
the same PARAMS values as dyn_rule 'slabs'. The only Earth inputs: trench points per Myr (trenches.py) and the
subducting floor's age there (Müller 2019 age grids, ages.npz). Steps t = 250 -> 0 Ma as the planet's: advance, deposit,
merge. Prints the air-loaded rms (l 1-31), its share in l 1-3, the deposition rate, the record by depth.
"""
import json
import sys
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import DATA
from planet_tectonics import mantle_flow as mf
from planet_tectonics.params import PARAMS

P = PARAMS
R = 6371e3


def main(trench_file, age_file, out_json, alpha='depth'):
    T = np.load(trench_file); A = np.load(age_file)
    age = A['age_myr']
    assert age.shape == T['t'].shape
    g = CubeSphere(217)
    z_cut = P['flow_z_cut_km']*1e3; lmax = int(P['flow_lmax'])
    radii = mf.load_radii(R, z_cut)
    base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08')
    K = mf.FlowKernels(R, radii, base, eta, lmax)
    factor = mf.alpha_ratio((R-radii)/1e3) if alpha == 'depth' else None
    rec = mf.SlabRecord(P['slab_dip_deg'], P['lm_sink_mm_yr'])
    out = {'series': [], 'inputs': {'trenches': str(trench_file), 'ages': str(age_file), 'alpha': alpha}}
    tt = T['t']
    for t in range(250, -1, -1):
        rec.advance(1., R)
        m = (tt == t) & np.isfinite(age) & (T['conv_mm_yr'] > 0)
        conv = T['conv_mm_yr'][m]*1e-3                                   # m/yr
        mass = mf.thermal_mass_excess(age[m], P['rho_mantle'], P['rho_water'])*conv*1e6*T['length_km'][m]*1e3
        rec.deposit(T['xyz'][m], T['n_over'][m], mass, conv, -float(t))
        rec.merge(g, z_cut)
        if t % 10 == 0 or t <= 20:
            loads = rec.layer_loads(z_cut, len(radii), R, lmax, 1, factor)
            H, V = K.surface(loads)
            Pl = mf.spectrum(H/P['rho_mantle'], lmax)
            row = {'t_ma': t, 'rms_m': float(np.sqrt(Pl[1:].sum())), 'share_l1_3': float(Pl[1:4].sum()/Pl[1:].sum()),
                   'deg_rms_m': {l: float(np.sqrt(Pl[l])) for l in (1, 2, 3, 4, 6, 8, 12, 20, 30)},
                   'flow_rms_cm_yr': float(np.sqrt(mf.gradient_power(V*3.15576e7, lmax))*100),
                   'deposit_kg_yr': float(mass.sum()/1e6), 'consumed_km2_yr': float((conv*T['length_km'][m]*1e3).sum()/1e6),
                   'held_below_cut_kg': float(rec.mass_by_layer(z_cut, len(radii)).sum())}
            out['series'].append(row)
            print(f"t={t:4d} Ma: rms {row['rms_m']:.0f} m (l1-3 {row['share_l1_3']:.2f}), flow {row['flow_rms_cm_yr']:.2f} cm/yr, "
                  f"deposit {row['deposit_kg_yr']:.2e} kg/yr, consumed {row['consumed_km2_yr']:.2f} km2/yr", flush=True)
    mbl = rec.mass_by_layer(z_cut, len(radii))
    zc = z_cut/1e3+25+50*np.arange(len(mbl))
    bins = [(300, 660), (660, 1000), (1000, 1500), (1500, 2000), (2000, 2500), (2500, 2850)]
    out['record_by_bin_1e21_kg_surface_alpha'] = {f'{a}-{b}': float(mbl[(zc > a) & (zc < b)].sum()/1e21) for a, b in bins}
    out['record_ledger_kg'] = {'deposited': rec.deposited_kg, 'held': rec.total_mass(), 'exited': rec.exited_kg}
    np.savez_compressed(out_json.replace('.json', '_H.npz'), H=H/P['rho_mantle'], V=V, lmax=lmax)
    json.dump(out, open(out_json, 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ('record_by_bin_1e21_kg_surface_alpha', 'record_ledger_kg')}, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:4], *(sys.argv[4:5] or []))
