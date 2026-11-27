"""Pass 59 probe: where the plate layer's export comes from, against Earth under the same law (journal 81).

One plate step (dt_myr 1, erosion_substeps 4) of the plate layer's own 46 km erosion (PlateSim._erode's loop:
graph_landscape.surface_step + fill_and_spill + shelf_step, PARAMS' law: K x 1.0, theta 200, Sc 0.036, Kd x 4) on
  Earth: ETOPO on n = 217 (outputs/planet-reference/earth-n217.npz); erodibility 1 on Bird 2003's orogens, 0.15
         elsewhere (as scripts.erosion_law_check); sediment 0 (only the shelf step reads it);
  planet: each run's present snapshot (t+0000.0.npz): elevation_m (relative to sea level), continental,
         orogeny_age_myr (erodibility 0.15 + 0.85 exp(-age / 50 Myr) on continental nodes, 1 on oceanic), sediment_m.
Measures, land = nodes above sea level at the start of the step, area-weighted:
  G gross erosion on land (km^3/yr); X net loss of land volume (eroded - deposited on land: Earth's sense of river
  delivery, as the runs count export); by class of 46 km slope (|grad z| on the graph, m/km), of erodibility (young
  orogen >= 0.5 / stable) and of elevation.
Shift-share on the slope x erodibility classes c (planet p, Earth E; r = G per area):
  G_p - G_E = sum_c (A_p,c - A_E,c) r_E,c   [terrain: how much land of each kind]
            + sum_c A_E,c (r_p,c - r_E,c)   [rate at equal slope and erodibility]
            + sum_c (A_p,c - A_E,c)(r_p,c - r_E,c)   [interaction]
  and the terrain term split again into slope x elevation-class areas at fixed erodibility.
Rule (fixed 2 Nov before the probe was run; journal 81 section 3), six-seed means of combo-02 (the export base:
C1 failed on combo-03), combo-04 and the control (sweep-20261101-115244 combo-01) reported beside:
  sanity: Earth's X under the law inside 4-8 km^3/yr (ice-sheet regions out) and the planet's one-step X within
          +-35% of its run's export (elevation_check); else the probe is not read.
  the term (terrain or rate) carrying >= 50% of the excess G_p - G_E is the export's cause;
    terrain -> the next step is the belts' shape (how much land of which slope and erodibility: item 1's crust
               distribution), not a rate;
    rate    -> the law's inputs at equal slope and erodibility differ from Earth's: the erodibility clock (orogen_rule
               belt resets continually on thick crust) is next, read against Earth's orogens.
  Within the terrain term, the class pair that carries the most is named (steep young belts or broad high land).
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_tectonics import shelf as shelf_mod
from planet_tectonics.params import PARAMS as P

SLOPE_EDGES = (0., 2., 5., 10., 20., 40., 1e9)          # m/km
ELEV_EDGES = (0., 200., 500., 1000., 2000., 3000., 1e9)  # m


def law_step(g, z0, cont, erod, sed):
    """PlateSim._erode's loop on surface z0 (relative to sea level). Returns dz (m over one plate step)."""
    p = P
    lp = {**gl.DEFAULTS, 'K': gl.DEFAULTS['K']*p['erosion_K_factor'], 'Sc': p['erosion_Sc'],
          'Kd': gl.DEFAULTS['Kd']*p['erosion_Kd_factor'], 'open_edges': False, 'theta': p['erosion_theta']}
    K = lp['K']*np.where(cont, erod, 1.)
    seed = z0 < -1500.
    dt = p['dt_myr']*1e6/p['erosion_substeps']
    z = z0.astype(np.float64).copy()
    fb = np.ones(g.size, bool)
    h_min, bound = g.spacing_min(), g.diffusion_bound()
    for _ in range(p['erosion_substeps']):
        r = gl.surface_step(g, z, seed, None, fb, K, lp, dt, h_min, bound, keep_deep=p['keep_deep_sediment'])
        z = r['z']
        if p['keep_deep_sediment'] and r['deep_m3'].any():
            z = gl.fill_and_spill(z, r['deep_m3'], g.nbr, g.area)
        if p['shelf_kappa_m2_yr'] > 0:
            ocean = gl.connected_ocean(z, seed, g.nbr)
            avail = np.where(cont, 1e9, np.clip(sed+z-z0, 0, None))
            z, _, _ = shelf_mod.shelf_step(g, z, ocean, avail, p['shelf_kappa_m2_yr'], p['shelf_wave_base_m'],
                                           p['cliff_retreat_m_yr'], dt)
    return z-z0


def slope_m_km(g, z):
    """Steepest descent to a neighbour, m/km (the law's own measure of the local gradient)."""
    nb = np.clip(g.nbr, 0, None)
    d = (z[:, None]-z[nb])/(g.dist/1e3)
    d = np.where(g.nbr >= 0, d, 0.)
    return np.clip(d.max(1), 0, None)


def classes(g, z0, erod, dz, land):
    a = g.area
    s = slope_m_km(g, z0)
    si = np.digitize(s, SLOPE_EDGES[1:-1])
    ei = np.digitize(np.clip(z0, 0, None), ELEV_EDGES[1:-1])
    oi = (erod >= .5).astype(int)
    eroded = np.clip(-dz, 0, None)
    out = {'G_km3_yr': float((eroded*a)[land].sum()/1e9/1e6), 'X_km3_yr': float((-dz*a)[land].sum()/1e9/1e6),
           'land_area_mkm2': float(a[land].sum()/1e12), 'cells': {}}
    for o in (0, 1):
        for i in range(len(SLOPE_EDGES)-1):
            for j in range(len(ELEV_EDGES)-1):
                m = land & (oi == o) & (si == i) & (ei == j)
                out['cells'][f'{o}{i}{j}'] = (float(a[m].sum()/1e12), float((eroded*a)[m].sum()/1e9/1e6))
    return out


def earth(g, ice_out=True):
    z = np.load(Path('outputs/planet-reference/earth-n217.npz'))['elevation_m'].astype(np.float64)
    sys.path.insert(0, '.')
    from planet_earth.bird import orogen_mask
    oro = orogen_mask(g, 'data/tectonicplates-master/GeoJSON/PB2002_orogens.json') if True else None
    erod = np.where(oro, 1., P['stable_erodibility'])
    cont = z > -1500.
    dz = law_step(g, z, cont, erod, np.zeros(g.size))
    land = z > 0
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1)))
    lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ice = (lat < -60) | ((lat > 59) & (lon > -75) & (lon < -10) & (lat < 84))
    res = {'all': classes(g, z, erod, dz, land)}
    res['no_ice'] = classes(g, z, erod, dz, land & ~ice)
    return res


def planet(g, f):
    s = np.load(f)
    z = s['elevation_m'].astype(np.float64)
    cont = s['continental'].astype(bool)
    age = s['orogeny_age_myr'].astype(np.float64)
    since = np.where(np.isnan(age), np.inf, age)
    erod = P['stable_erodibility']+(1-P['stable_erodibility'])*np.exp(-np.clip(since, 0, None)/P['orogen_erodibility_tau_myr'])
    erod_k = np.where(cont, erod, 1.)
    dz = law_step(g, z, cont, erod, s['sediment_m'].astype(np.float64))
    return classes(g, z, erod_k, dz, z > 0)


if __name__ == '__main__':
    g = CubeSphere(217)
    out = {'earth': earth(g)}
    for f in sys.argv[1:]:
        key = '/'.join(Path(f).parts[-5:-2])
        out[key] = planet(g, f)
        print(key, round(out[key]['G_km3_yr'], 2), round(out[key]['X_km3_yr'], 2), flush=True)
    Path('export_source.json').write_text(json.dumps(out))
    print('earth', {k: (round(v['G_km3_yr'], 2), round(v['X_km3_yr'], 2)) for k, v in out['earth'].items()})
