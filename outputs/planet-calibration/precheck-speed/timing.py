"""Seconds per step of a plate run as the driver steps it (step + state), journal 154.
  python timing.py n seed threads steps   (run from the checkout whose code is timed; PLANET_DATA for data/)"""
import sys, glob, time
sys.path[:0] = glob.glob('packages/*/src')
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim
n, seed, T, K = (int(v) for v in sys.argv[1:5])
g = CubeSphere(n)
sim = PlateSim(g, seed=seed, t0=-250., params={'threads': T or None})
sim.step(); sim.state()                                   # first step (one-off tables) not timed
t0 = time.perf_counter()
for _ in range(K):
    sim.step(); sim.state()
print(f'{(time.perf_counter()-t0)/K:.3f} s/step over {K} steps, seed {seed}')
