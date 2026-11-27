"""Pass 59, descriptive: land the law treats as young orogen (planet: continental, orogeny < 50 Myr; Earth: Bird's
orogens), its area and 46 km slope, six seeds per arm against Earth on the same grid. Reported, not scored."""
import glob, sys, json
import numpy as np
sys.path.insert(0, '/home/claude/probe')
from export_source import slope_m_km
from planet_earth.bird import orogen_mask
from planet_grid.sphere_graph import CubeSphere
g = CubeSphere(217); a = g.area
U = '/mnt/user-data/uploads/Planetary_Diffusion/outputs/planet-calibration/'
def desc(z, belt, land):
    s = slope_m_km(g, z)
    b = belt & land; o = land & ~belt
    return dict(belt_mkm2=a[b].sum()/1e12, other_mkm2=a[o].sum()/1e12, belt_slope_p50=float(np.median(s[b])), belt_slope_p90=float(np.percentile(s[b], 90)),
                other_slope_p50=float(np.median(s[o])), belt_z_p50=float(np.median(z[b])), other_z_p50=float(np.median(z[o])))
z = np.load('outputs/planet-reference/earth-n217.npz')['elevation_m'].astype(float)
oro = orogen_mask(g, 'data/tectonicplates-master/GeoJSON/PB2002_orogens.json')
lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
ice = (lat < -60) | ((lat > 59) & (lon > -75) & (lon < -10) & (lat < 84))
print('Earth', {k: round(v, 2) for k, v in desc(z, oro, (z > 0) & ~ice).items()})
for arm in ('sweep-20261101-191306/combo-02', 'sweep-20261101-191306/combo-03', 'sweep-20261101-191306/combo-04', 'sweep-20261101-115244/combo-01'):
    rows = []
    for f in sorted(glob.glob(U+arm+'/seed*/snapshots/t+0000.0.npz')):
        s = np.load(f); z = s['elevation_m'].astype(float); age = s['orogeny_age_myr'].astype(float)
        belt = s['continental'] & (np.nan_to_num(age, nan=1e9) < 50)
        rows.append(desc(z, belt, z > 0))
    print(arm, {k: round(float(np.mean([r[k] for r in rows])), 2) for k in rows[0]})
