import numpy as np, time
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import plates as P
g = CubeSphere(217); a = g.area
E = np.load('/home/claude/pre/earth_nodes.npz'); T = E['T']; cont = T >= 20e3
print('Earth vol %.0f Mkm3' % ((T*a)[cont].sum()/1e15))
bins = [0, 50, 100, 150, 200, 300, 400, 600, 800, 1000, 1500, 2000, 3000, 1e9]
t0 = time.time()
for seed in (1011, 1001):
    sim = P.PlateSim(g, seed=seed, params={'erosion': False})
    c = sim.world['cont']; th = sim.world['thick'].astype(float)
    dp = gl.distance_to(~c, g.nbr, g.dist, 2e7)/1e3
    X = th
    print('seed', seed, 'init %.0fs' % (time.time()-t0), 'area %.1f vol %.0f' % (a[c].sum()/1e12, (th*a)[c].sum()/1e15), 'dmax %.0f' % dp[c].max(), 'X20-30 %.1f' % a[c & (X < 30e3)].sum()/1e12 if False else '')
    print('   X20-30 %.1f' % (a[c & (X < 30e3)].sum()/1e12))
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = c & (dp > lo) & (dp <= hi)
        if m.any(): print('   %5d-%5d %6.1f  p50 %.1f' % (lo, min(hi, 9999), a[m].sum()/1e12, np.median(X[m])/1e3))
