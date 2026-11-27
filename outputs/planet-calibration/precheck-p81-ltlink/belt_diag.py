"""Pass 81: why 'shorten' leaves young thick lithosphere thin. n = 96, seed 1012, belt_lithosphere 'shorten', 120 steps:
per step after 60, the belt's continental nodes, those margin shortening thickened (belt_gain > 0), and young (< 250 Myr)
lithosphere under crust > 45 km split by what touched it this step: shortened, reset by the belt, neither."""
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.lithosphere import thickness_km
g = CubeSphere(96)
sim = P.PlateSim(g, seed=1012, params={'threads': 1, 'belt_lithosphere': 'shorten', 'start_lithosphere': 'crust'})
sim.state()
acc = np.zeros(6)
for k in range(120):
    sim.step(); st = sim.state()
    if k < 60:
        continue
    w = sim.world; c = w['cont']
    belt = c & sim.belt; gain = getattr(sim, 'belt_gain', np.zeros(g.size)) > 0
    T = st['crust_thickness_m']; la = st['lithosphere_age_myr']
    yt = c & (T > 45e3) & (la < 250)
    acc += [belt.sum(), (belt & gain).sum(), yt.sum(), (yt & gain).sum(), (yt & belt & ~gain).sum(), (yt & ~belt).sum()]
    if k % 15 == 0:
        print(f'step {k}: belt {belt.sum()}, shortened {(belt & gain).sum()}; young thick {yt.sum()} (shortened {(yt & gain).sum()}, '
              f'belt reset {(yt & belt & ~gain).sum()}, off the belt {(yt & ~belt).sum()}); L young thick '
              f'{thickness_km(la[yt], sim.p).mean() if yt.any() else float("nan"):.0f} km; young-thick age p50 '
              f'{np.median(la[yt]) if yt.any() else float("nan"):.1f} Myr')
b, bg, yt, ytg, ytr, yto = acc
print(f'steps 60-119: belt nodes shortened {bg/b:.2f}; young thick nodes: shortened {ytg/yt:.2f}, reset by the belt {ytr/yt:.2f}, off the belt {yto/yt:.2f}')
