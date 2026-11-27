import numpy as np, time
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_earth.crust1 import C1, crust1_cell
from planet_reviews.crust_review import wpct
g = CubeSphere(217); a = g.area; cell = crust1_cell(g)
T = np.loadtxt(C1/'crsthk.xyz')[:, 2][cell]*1e3; S = np.loadtxt(C1/'sedthk.xyz')[:, 2][cell]*1e3
cont = T >= 20e3; X = T-S
d = gl.distance_to(~cont, g.nbr, g.dist, 1e7)/1e3
print('Earth label-free area %.1f Mkm2 frac %.4f vol %.0f Mkm3 sed %.1f Mkm3' % (a[cont].sum()/1e12, a[cont].sum()/a.sum(), (T*a)[cont].sum()/1e15, (S*a)[cont].sum()/1e15))
print('d max', d[cont].max())
bins = [0, 50, 100, 150, 200, 300, 400, 600, 800, 1000, 1500, 2000, 3000, 1e9]
print('d band   area   X p10 p50 p90   S mean  X20-30 area')
for lo, hi in zip(bins[:-1], bins[1:]):
    m = cont & (d > lo) & (d <= hi)
    if not m.any(): continue
    print('%5d-%5d %6.1f  %5.1f %5.1f %5.1f  %6.0f  %5.1f' % (lo, min(hi, 9999), a[m].sum()/1e12, wpct(X[m], a[m], 10)/1e3, wpct(X[m], a[m], 50)/1e3, wpct(X[m], a[m], 90)/1e3, (S[m]*a[m]).sum()/a[m].sum(), a[m & (X >= 20e3) & (X < 30e3)].sum()/1e12))
np.savez('/home/claude/pre/earth_nodes.npz', T=T, S=S, d=d)
# the planet's start by distance: planet_by_d.py
