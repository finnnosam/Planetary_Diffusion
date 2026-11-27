"""Pass 52, before the rule: the measurements alone (no model number). Holdt et al. 2022 oceanic 1 deg bins."""
import numpy as np
from pathlib import Path
from planet_grid.paths import DATA
D = DATA/'residual_topography'/'Residual depthtopographic measurements'/'Residual depthtopographic measurements'/'Holdt_et_al_2022_SupportingInformation'
def load(n):
    return np.loadtxt(D/f'oceanic_measurements_{n}_1_deg.xyz')
for n in ('circles', 'down_triangles', 'up_triangles'):
    a = load(n); v = a[:, 2]; e = 0.5*(a[:, 4]-a[:, 3])
    blk = np.unique(np.floor((a[:, 1]+90)/10)*36+np.floor(np.mod(a[:, 0], 360)/10))
    print(f'{n:15s} n {len(v):5d} mean {v.mean():+.3f} km, std {v.std():.3f}, rms {np.sqrt((v**2).mean()):.3f}, '
          f'p5/p50/p95 {np.percentile(v, 5):+.2f}/{np.median(v):+.2f}/{np.percentile(v, 95):+.2f}, mean err {e.mean():.3f}, '
          f'10-deg blocks {blk.size}, lat range {a[:,1].min():.0f}..{a[:,1].max():.0f}')
c = np.loadtxt(D/'continental_measurements_seismic_1_deg.xyz')[:, 2]
print(f'continental seismic n {c.size} mean {c.mean():+.3f} std {c.std():.3f}')
