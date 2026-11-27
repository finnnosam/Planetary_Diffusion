"""Pass 50, reported only: the load-mass scale s that the geoid prefers (least squares, s = sum(obs.mod)/sum(mod^2)) per
degree band, against the scale the observed residual topography needs (220 m / E3's rms), per candidate profile."""
import json
import numpy as np
import profiles as p
from planet_tectonics import mantle_flow as mf
out = {}
for name, (base, eta) in p.PROFILES.items():
    K = mf.FlowKernels(p.A, p.RAD, base, eta, p.LMAX)
    loads = p.E3['loads_0']*p.FAC[:, None]
    H, _ = K.surface(loads); N = K.geoid(loads)
    topo = float(np.sqrt(mf.spectrum(H/p.RHO_M, p.LMAX)[1:].sum()))
    row = {'s_topo_0.22': 220./topo, 's_topo_0.15_0.30': [150./topo, 300./topo]}
    for band in ((2, 8), (4, 9), (2, 3)):
        s = np.isin(p.li, range(band[0], band[1]+1)); w = p.W[s]
        for k, o in p.OBS.items():
            sc = float((w*(o[s]*np.conj(N[s])).real).sum()/(w*np.abs(N[s])**2).sum())
            row[f's_geoid_l{band[0]}_{band[1]}_{k}'] = sc
            row[f'VR_at_s_topo_l{band[0]}_{band[1]}_{k}'] = p.vr(N*220./topo, o, range(band[0], band[1]+1))
            row[f'VR_l{band[0]}_{band[1]}_{k}'] = p.vr(N, o, range(band[0], band[1]+1))
    out[name] = row
    print(f"{name}: topo needs s = {row['s_topo_0.22']:.2f} ({row['s_topo_0.15_0.30'][0]:.2f}-{row['s_topo_0.15_0.30'][1]:.2f}); geoid prefers "
          f"l4-9 {row['s_geoid_l4_9_chambat']:.2f}, l2-8 {row['s_geoid_l2_8_chambat']:.2f}/{row['s_geoid_l2_8_nakiboglu']:.2f}, l2-3 {row['s_geoid_l2_3_chambat']:.2f}; "
          f"VR l4-9 at s=1 {row['VR_l4_9_chambat']:+.2f}, at the topography's s {row['VR_at_s_topo_l4_9_chambat']:+.2f}; "
          f"VR l2-8 at the topography's s {row['VR_at_s_topo_l2_8_chambat']:+.2f}/{row['VR_at_s_topo_l2_8_nakiboglu']:+.2f}")
json.dump(out, open('massscale.json', 'w'), indent=1)
