"""Pass 103 probe: the flow's traction under the plate holding the most continent, at 0 Myr, from the snapshots (the
velocity the rift_traction arm would read: flow_velocity_alm_m_yr, each degree x sqrt(target / raw power), target the
measured power x flow_placed's share). Per seed: the plate's continent, rms flow speed under it (cm/yr), the mean
outward speed from its centroid (cm/yr; > 0 divergent), and the traction's outward force summed over the plate
(C_b (1 + (keel-1) c) a v.r_hat, N; r_hat away from the centroid) beside one 1,000 km trench's slab pull (F_sp x 1e6 m). Read. Usage: python3 traction_read.py <arm>"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_tectonics.mantle import measured_power_air_m2, placed_share
from planet_tectonics.plates import PARAMS
from planet_reviews.crust_review import sediment_bulk_m
g = CubeSphere(217); a = g.area
arm = Path(sys.argv[1])
for s in sorted(int(p.name[4:]) for p in arm.glob("seed*")):
    q = np.load(arm/f'seed{s}'/'snapshots'/'t+0000.0.npz')
    lm = int(q['flow_lmax']); L = mf.alm_l(lm)
    tp = np.zeros(lm+1); pw = measured_power_air_m2(); tp[:min(pw.size, lm+1)] = pw[:lm+1]; tp[0] = 0
    tp = tp*placed_share(lm)
    raw = q['flow_power_raw_m2']; sc = np.where(raw > 0, np.sqrt(tp/np.where(raw > 0, raw, 1)), 0)
    v = mf.gradient(q['flow_velocity_alm_m_yr'].astype(np.complex128)*sc[L], lm, g.xyz)
    c = (q['crust_thickness_m']+sediment_bulk_m(q)) >= 20e3
    pid = q['plate_id']; ids, inv = np.unique(pid, return_inverse=True)
    pc = np.bincount(inv, weights=a*c); k = np.argmax(pc); own = inv == k
    cen = (a[own, None]*g.xyz[own]).sum(0); cen /= np.linalg.norm(cen)
    r = cen[None, :]-(g.xyz@cen)[:, None]*g.xyz; r = -r/np.maximum(np.linalg.norm(r, axis=1), 1e-12)[:, None]
    out = (v*r).sum(1)
    w = a*(1+(PARAMS['keel_drag_factor']-1)*c)
    F = PARAMS['basal_drag_pa_s_m']*(w[own]*out[own]).sum()/3.15576e7
    print(f'seed {s}: plate continent {pc[k]/1e12:5.0f} Mkm2, area {a[own].sum()/1e12:5.0f}; flow rms under it '
          f'{100*np.sqrt(np.average((v[own]**2).sum(1), weights=a[own])):.2f} cm/yr, mean outward {100*np.average(out[own], weights=a[own]):+.2f} cm/yr; '
          f'outward traction {F:.2e} N; one 1,000 km trench pull {PARAMS["slab_pull_n_m"]*1e6:.2e} N; '
          f'raw-to-placed factor l 1/2/3 {sc[1]:.2f}/{sc[2]:.2f}/{sc[3]:.2f}')
