"""Pass 101: E3 (Earth's subduction history through D, journal 70) with today's PARAMS: the present-day load layers
saved in pass 50 (precheck-p50-geoid/e3_loads.npz, no expansivity factor) times slab_alpha 'thermo' (mantle_flow.
thermo_ratio, the default since pass 51), through MantleFlow's kernels (visc.sh08, cut 300 km, l 1-31). Writes
e3_thermo_H.npz (H air-loaded m, lmax) and prints its rms and the 'depth' version's for comparison.
Usage (repo root): python3 e3_thermo.py"""
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_tectonics import mantle_flow as mf
from planet_tectonics.params import PARAMS as P
d = np.load(OUTPUTS/'planet-calibration'/'precheck-p50-geoid'/'e3_loads.npz')
R = 6371e3; radii = d['radii']; lmax = int(d['lmax'])
base, eta = mf.read_viscosity(DATA/'mantle'/'visc.sh08')
K = mf.FlowKernels(R, radii, base, eta, lmax)
out = {}
for name, fac in (('thermo', mf.thermo_ratio((R-radii)/1e3)), ('depth', mf.alpha_ratio((R-radii)/1e3))):
    H, V = K.surface(d['loads_0']*fac[:, None])
    H = H/P['rho_mantle']; out[name] = H
    print(f'{name}: rms {np.sqrt(mf.spectrum(H, lmax)[1:].sum()):.0f} m')
H0 = np.load(OUTPUTS/'planet-calibration'/'precheck-p49-e3'/'e3_H.npz')['H']
print(f'saved e3_H rms {np.sqrt(mf.spectrum(H0, lmax)[1:].sum()):.0f} m; max |H_depth - e3_H| {np.abs(out["depth"]-H0).max():.2e} m')
np.savez(Path(__file__).with_name('e3_thermo_H.npz'), H=out['thermo'], lmax=lmax)
