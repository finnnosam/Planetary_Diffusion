import numpy as np, glob, re
from planet_grid.sphere_graph import CubeSphere
g = CubeSphere(217); A = g.area/g.area.sum()
fs = sorted(glob.glob('outputs/planet-plates/freeze-20261031/seed1001/snapshots/t*.npz'), key=lambda f: float(re.search(r't([+-][\d.]+)\.npz', f).group(1)))
print(' t_myr  cont  >45/cont  >55/cont  T_cont_p50  vol_rel  col_share  col_T_p50  >45 by orogeny age <25/25-75/75-150/>150 Myr (share of cont)')
v0 = None
def wm(v, w): o = np.argsort(v); c = np.cumsum(w[o]); return float(v[o][np.searchsorted(c, c[-1]/2)])
for f in fs:
    t = float(re.search(r't([+-][\d.]+)\.npz', f).group(1))
    if t > -40 and t % 10: continue
    s = np.load(f); T = s['crust_thickness_m'].astype(float); sed = s['sediment_m'].astype(float); tot = T+sed
    c = tot >= 20e3; ca = A[c].sum(); vol = (A*T*s['continental']).sum()
    v0 = vol if v0 is None else v0
    oa = s['orogeny_age_myr'].astype(float); runlen = t+250
    col = c & (oa < runlen)
    thick = c & (T > 45e3)
    bands = [thick & (oa < 25), thick & (oa >= 25) & (oa < 75), thick & (oa >= 75) & (oa < 150), thick & (oa >= 150)]
    print(f'{t:6.0f} {ca:.3f}  {A[thick].sum()/ca:.3f}    {A[c & (T>55e3)].sum()/ca:.3f}     {wm(T[c],A[c])/1e3:5.1f}     {vol/v0:.3f}   {A[col].sum()/ca:.2f}      {wm(T[col],A[col])/1e3 if col.any() else float("nan"):5.1f}     ' + ' / '.join(f'{A[b].sum()/ca:.3f}' for b in bands))
