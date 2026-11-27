"""Pass 83 (journal 108): S's room on 0 Myr states. Where S applies (off the label-free continent, > 300 km from any
plate boundary, within 1,000 km of the continent), the volume between the planet's floor and Earth's measured rise
profile (PARAMS margin_wedge_profile), as elevation (Mkm^3) and as the sediment that fills it (x (rho_m - rho_w) /
(rho_m - rho_s), bulk 2,000 kg/m^3), and the floor's mean offset from the profile by distance.
Usage (repo root): python wedge_room.py <arm> [<arm> ...]"""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics.params import PARAMS
g = CubeSphere(217); a = g.area; AREA = 4*np.pi*6371e3**2
km, depth = np.array(PARAMS['margin_wedge_profile']).T
BINS = ((0, 100), (100, 200), (200, 300), (300, 500), (500, 700), (700, 1000))
for arm in sys.argv[1:]:
    rooms, offs = [], []
    for s in range(1001, 1007):
        d = np.load(f'{arm}/seed{s}/snapshots/t+0000.0.npz')
        z = d['elevation_m'].astype(float); col = d['crust_thickness_m'].astype(float)+d['sediment_bulk_m'].astype(float)
        edge = col >= 20e3
        dc = gl.distance_to(edge, g.nbr, g.dist, 1e6)
        far = gl.distance_to(d['boundary_type'] != 0, g.nbr, g.dist, 3e5) > 3e5
        on = ~edge & far & (dc <= 1e6)
        top = np.interp(dc/1e3, km, depth)
        r = np.where(on, np.clip(top-z, 0, None), 0.)
        rooms.append((r*a).sum()*AREA/a.sum()/1e9/1e6)
        offs.append([((z-top)[on & (dc/1e3 >= lo) & (dc/1e3 < hi)]*a[on & (dc/1e3 >= lo) & (dc/1e3 < hi)]).sum()
                     / max(a[on & (dc/1e3 >= lo) & (dc/1e3 < hi)].sum(), 1) for lo, hi in BINS]
                    + [a[on].sum()/a.sum()])
    rooms = np.array(rooms); offs = np.array(offs)
    print(Path(arm).name, f': room {rooms.mean():.1f} Mkm3 of elevation (seeds {np.round(rooms, 1)}), '
          f'~{rooms.mean()*2270/1300:.0f} Mkm3 of bulk sediment to fill it; S applies on {offs[:, -1].mean():.3f} of the surface')
    print('  floor minus profile (m) by distance', ' '.join(f'{lo}-{hi}: {v:+.0f}' for (lo, hi), v in zip(BINS, offs[:, :-1].mean(0))))
