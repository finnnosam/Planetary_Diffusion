"""Pass 59 probe: where the plate layer erodes, against Earth's measured erosion (journal 81 section 4).

Planet: erosion_record.npz from erosion_record_run.py (cloud seed 1009, a seventh read: not pooled, not scored on the
six seeds). Over the last 40 Myr, per node: gross continental erosion E = -min(e_cont + dsed, 0) summed per step
(crust and continental sediment removed, m), this step's tectonic thickening, the surface the law saw.
Classes (land: z0 > 0 in at least half of the recorded steps):
  belt   continental nodes with an orogeny younger than 50 Myr in at least half the steps (erodibility >= 0.46: what
         the law treats as a young orogen);
  other  the rest of the land.
Earth: Bird 2003's orogens (PB2002_orogens) on ETOPO land on n = 217, ice-sheet regions out; rates from OCTOPUS v2
basins at least one 46 km cell in area (710 basins, outputs/erosion-law-check/217-20261023-202816/results.json:
measured rate, area, centre), each class's rate the area-weighted mean of its basins (centre node's class).
Rule (fixed 2 Nov, before either record was read):
  G = sum_c A_c r_c (km^3/yr), planet against Earth; the excess split
    area term  sum_c (A_p,c - A_E,c) r_E,c;   rate term  sum_c A_E,c (r_p,c - r_E,c);   interaction the rest.
  The term carrying >= 50% of the excess names the cause of the export:
    area -> too much land the law treats as young orogen (the orogeny clock / the crust distribution);
    rate -> the belts (or the rest) erode faster than Earth's at the same class: their steepness at 46 km.
  Planet E (amended 2 Nov before reading: the record sums per-step changes) = crystalline erosion + net sediment
  lost over the 40 Myr: a lower bound on the gross, so a rate excess read from it is not overstated.
  Supply beside it: positive tectonic thickening on belts against erosion there (steady state: about equal).
  Both arms are read the same way (base: combo-02's settings; params: today's PARAMS + measured field).
"""
import json
import sys
from pathlib import Path

import numpy as np

from planet_earth.bird import orogen_mask
from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere

TAU = 50.


def earth(g):
    z = np.load('outputs/planet-reference/earth-n217.npz')['elevation_m'].astype(np.float64)
    oro = orogen_mask(g, 'data/tectonicplates-master/GeoJSON/PB2002_orogens.json')
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1)))
    lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    ice = (lat < -60) | ((lat > 59) & (lon > -75) & (lon < -10) & (lat < 84))
    land = (z > 0) & ~ice
    rows = json.loads(Path(sys.argv[1]).read_text())['rows']
    lo = np.array([r['lon'] for r in rows]); la = np.array([r['lat'] for r in rows])
    xyz = np.stack((np.cos(np.radians(la))*np.cos(np.radians(lo)), np.cos(np.radians(la))*np.sin(np.radians(lo)),
                    np.sin(np.radians(la))), -1)
    node = fields.nearest_nodes(g, xyz)
    rate = np.array([r['measured_m_myr'] for r in rows]); area = np.array([r['area_km2'] for r in rows])
    out = {}
    for name, m in (('belt', land & oro), ('other', land & ~oro)):
        b = m[node]
        out[name] = {'area_mkm2': float(g.area[m].sum()/1e12), 'basins': int(b.sum()),
                     'rate_m_myr': float((rate[b]*area[b]).sum()/area[b].sum()),
                     'rate_median_m_myr': float(np.median(rate[b]))}
    return out


def planet(g, f):
    d = np.load(f)
    n = int(d['steps'])
    a = d['area']
    # the record sums per-step changes: crust erosion (e_cont <= 0) is gross; sediment (dsed +-) only nets out, so
    # E = crust eroded + net sediment lost where the node lost sediment over the 40 Myr (a lower bound on the gross)
    E = -d['acc_e_cont']+np.clip(-d['acc_dsed'], 0, None)
    land = d['acc_land_steps'] >= n/2
    belt = land & (d['acc_orog_steps'] >= n/2)
    other = land & ~belt
    out = {'steps': n}
    for name, m in (('belt', belt), ('other', other)):
        out[name] = {'area_mkm2': float(a[m].sum()/1e12), 'rate_m_myr': float((E*a)[m].sum()/a[m].sum()/n),
                     'tect_pos_m_myr': float((d['acc_tect_pos']*a)[m].sum()/a[m].sum()/n),
                     'tect_net_m_myr': float((d['acc_tect']*a)[m].sum()/a[m].sum()/n)}
    out['submerged_cont_G_km3_yr'] = float((E*a)[~land].sum()/n/1e15)
    return out


def split(p, e):
    G = lambda x: sum(x[c]['area_mkm2']*1e12*x[c]['rate_m_myr']/1e6/1e9 for c in ('belt', 'other'))
    gp, ge = G(p), G(e)
    at = sum((p[c]['area_mkm2']-e[c]['area_mkm2'])*1e12*e[c]['rate_m_myr']/1e15 for c in ('belt', 'other'))
    rt = sum(e[c]['area_mkm2']*1e12*(p[c]['rate_m_myr']-e[c]['rate_m_myr'])/1e15 for c in ('belt', 'other'))
    return {'G_planet': gp, 'G_earth': ge, 'excess': gp-ge, 'area_term': at, 'rate_term': rt, 'interaction': gp-ge-at-rt}


if __name__ == '__main__':
    g = CubeSphere(217)
    e = earth(g)
    print('Earth', json.dumps(e))
    res = {'earth': e}
    for f in sys.argv[2:]:
        p = planet(g, f)
        res[f] = {'planet': p, 'split': split(p, e)}
        print(f, json.dumps(p), json.dumps(res[f]['split']))
    Path('erosion_where.json').write_text(json.dumps(res, indent=1))
