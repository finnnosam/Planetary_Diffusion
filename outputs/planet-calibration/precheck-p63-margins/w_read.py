# Pass 63: journal 84's W readings on runs at 0 Myr (W1 crust 27.5-32.5 km label-free, W2 margin profile at 1 degree,
# W3 width to 35 km p50 at 1 degree, W4 drowned continent per km of coast), with the shelf, the sea on the datum, the
# label-free area and the continental volume change.  python w_read.py <seed dirs...>
import sys, json, numpy as np
from pathlib import Path
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params
from planet_reviews.margin_width import measure, to_crust1_cells

g = CubeSphere(217); a = g.area; nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(a.mean()))
print('run | W1 27.5-32.5 Mkm2 (>= 24.2) | W2 profile 50-100/100-150/150-200/200-300 (29.9/31.6/33.0/35.0 +-1) | interior 500-1000 '
      '(37.1) | W3 width35 p50 (>= 200) | W4 drowned per km coast (>= 100) | shelf Mkm2 (23.1) | passive k km (1 deg) | thin per km | '
      'area | sea m | volume change')
for s in map(Path, sys.argv[1:]):
    z = np.load(s/'snapshots'/'t+0000.0.npz'); p = run_params(s)
    T = column_thickness_m(z, p); cont = label_free_continent(T)
    zz = z['elevation_m'].astype(float); land = zz > 0; sea = ~land
    coast = land & np.any(np.where(ok, sea[np.where(ok, nb, 0)], False), 1)
    w1 = a[cont & (T >= 27.5e3) & (T < 32.5e3)].sum()/1e12
    T1 = to_crust1_cells(g, T, z['continental'].astype(bool))[0]
    m = measure(g, T1, label_free_continent(T1), z['boundary_type'] != 0)
    pr = m['profile_T_km_p25_50_75']
    w4 = a[cont & sea].sum()/1e6/(coast.sum()*sp/1e3)
    shelf = a[sea & (zz > -200)].sum()/1e12
    st = json.load(open(s/'stats.json'))
    vol = st.get('continental_volume_change', float('nan'))
    print(f"{s.parent.name}/{s.name}: W1 {w1:.1f} | W2 " + ' / '.join(f"{pr[k][1]:.1f}" for k in ('50-100', '100-150', '150-200', '200-300'))
          + f" | {pr['500-1000'][1]:.1f} | W3 {m['width_35_km_p25_50_75'][1]:.0f} | W4 {w4:.0f} | shelf {shelf:.1f} | passive "
          f"{m['passive_edge_length_km']/1e3:.1f} | thin/km {m['thin_per_km']:.0f} | area {a[cont].sum()/a.sum():.3f} | sea "
          f"{float(z['sea_level_m']):.0f} | vol {vol if isinstance(vol, float) else vol}")
