"""Pass 71: Earth's label-free continental pieces under 10 Mkm^2 (CRUST1.0 1-degree, n = 217): where, how big, edge and
passive edge (PB2002 within 300 km), and CRUST1.0's own crust type share (continental flag)."""
import sys, numpy as np
sys.path.insert(0, 'outputs/planet-calibration/precheck-p71-window')
from piece_edge import *
T, fl = earth_nodes(g); cont = label_free_continent(T)
edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
_, bd = nearest_label(np.where(isb, np.arange(g.size), -1).astype(np.int64), nb, g.dist); pas = edge & (bd > 300e3)
lab = comps(cont); A = np.bincount(np.maximum(lab, 0), weights=np.where(cont, a, 0), minlength=g.size)/1e12
for L in sorted(set(lab[cont]), key=lambda L: -A[L]):
    s = lab == L
    if A[L] >= 10: print(f"big piece {A[L]:.1f} Mkm2 edge {(edge&s).sum()*sp/1e6:.1f}k passive {(pas&s).sum()*sp/1e6:.1f}k cont-flag {fl[s].mean():.2f}"); continue
    c = g.xyz[s].mean(0); c /= np.linalg.norm(c); lat = np.degrees(np.arcsin(c[2])); lon = np.degrees(np.arctan2(c[1], c[0]))
    print(f"{lat:6.1f} {lon:7.1f} {A[L]:5.2f} Mkm2 edge {(edge&s).sum()*sp/1e3:5.0f} km passive {(pas&s).sum()*sp/1e3:5.0f} km cont-flag {fl[s].mean():.2f} T p50 {np.median(T[s])/1e3:.1f} km")
