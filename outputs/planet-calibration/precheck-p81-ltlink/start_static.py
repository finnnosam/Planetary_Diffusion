"""Pass 81: the start's L by crust bin (n = 217) under start_lithosphere 'distance' and 'crust', seeds given; beside
Earth's seismic L by crust bin. Run with PYTHONPATH at the repo."""
import sys, numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.lithosphere import thickness_km, elevation_term
g = CubeSphere(217)
EARTH = [97, 98, 107, 143, 158, 157, 197]
BINS = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]
for seed in map(int, sys.argv[1:]):
    res = {}
    for rule in ('distance', 'crust'):
        sim = P.PlateSim(g, seed=seed, params={'threads': 1, 'start_lithosphere': rule})
        st = sim.state()
        T = (st['crust_thickness_m']+st['sediment_m'])/1e3; la = st['lithosphere_age_myr']
        m = st['continental'] & np.isfinite(la) & (T >= 20)
        L = thickness_km(la[m], sim.p); a = g.area[m]; TT = T[m]
        res[rule] = ([np.average(L[(TT >= lo) & (TT < hi)], weights=a[(TT >= lo) & (TT < hi)]) if ((TT >= lo) & (TT < hi)).any() else np.nan
                      for lo, hi in BINS], np.polyfit(TT, L, 1)[0], np.corrcoef(TT, L)[0, 1], np.average(L, weights=a),
                     np.average(elevation_term(la[m], sim.p), weights=a))
    for rule, (b, s, r, mu, dz) in res.items():
        print(f'seed {seed} {rule:8s}: L by crust ' + ' / '.join(f'{x:.0f}' for x in b) + f'; dL/dT {s:.2f} corr {r:.2f}; L mean {mu:.1f}, C term mean {dz:.0f} m')
print('Earth seismic       : ' + ' / '.join(map(str, EARTH)) + ' (20-25 ... 55-80 km); dL/dT 3.18 on nodes')
