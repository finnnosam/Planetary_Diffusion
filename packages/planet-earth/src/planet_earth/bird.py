"""Readers for Bird (2003) PB2002 (data/tectonicplates-master): orogen polygons, plate-boundary lines and the
boundary-step table, on the cube-sphere graph."""
import json
import struct
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.paths import DATA


def lonlat(xyz):
    x, y, z = np.asarray(xyz).T
    return np.degrees(np.arctan2(y, x)), np.degrees(np.arcsin(np.clip(z, -1, 1)))


def orogen_mask(g, path):
    """Nodes inside Bird 2003's orogen polygons (longitudes unwrapped around each ring's first vertex)."""
    from matplotlib.path import Path as MPath
    lon, lat = lonlat(g.xyz)
    inside = np.zeros(g.size, bool)
    for f in json.loads(Path(path).read_text(encoding='utf-8'))['features']:
        geom = f['geometry']
        for poly in (geom['coordinates'] if geom['type'] == 'MultiPolygon' else [geom['coordinates']]):
            ring = np.asarray(poly[0], float)
            l0 = ring[0, 0]
            rl = (ring[:, 0]-l0+180) % 360-180+l0
            pl = (lon-l0+180) % 360-180+l0
            inside |= MPath(np.c_[rl, ring[:, 1]]).contains_points(np.c_[pl, lat])
    return inside


def earth_boundaries(g, step_km=10.):
    js = json.loads((DATA/'tectonicplates-master'/'GeoJSON'/'PB2002_boundaries.json').read_text())
    pts = []
    for f in js['features']:
        geo = f['geometry']
        lines = geo['coordinates'] if geo['type'] == 'MultiLineString' else [geo['coordinates']]
        for ln in lines:
            ll = np.radians(np.array(ln, float))
            xyz = np.stack((np.cos(ll[:, 1])*np.cos(ll[:, 0]), np.cos(ll[:, 1])*np.sin(ll[:, 0]), np.sin(ll[:, 1])), -1)
            for a, b in zip(xyz[:-1], xyz[1:]):
                n = max(2, int(np.ceil(np.arccos(np.clip(a@b, -1, 1))*g.radius/1e3/step_km))+1)
                s = np.linspace(0, 1, n)[:, None]
                p = a*(1-s)+b*s
                pts.append(p/np.linalg.norm(p, axis=1, keepdims=True))
    m = np.zeros(g.size, bool)
    m[fields.nearest_nodes(g, np.concatenate(pts))] = True
    return m


def read_dbf(path):
    """dBase III table (Bird's .dbf files) as {field: array}; numbers as float, text stripped."""
    raw = Path(path).read_bytes()
    n, hl, rl = struct.unpack('<IHH', raw[4:12])
    fl, o = [], 32
    while raw[o] != 0x0d:
        d = raw[o:o+32]
        fl.append((d[:11].split(b'\0')[0].decode(), chr(d[11]), d[16])); o += 32
    out = {k: [] for k, _, _ in fl}
    for i in range(n):
        r, q = raw[hl+i*rl+1:hl+(i+1)*rl], 0
        for k, t, w in fl:
            s = r[q:q+w].decode('latin1').strip(); q += w
            out[k].append((float(s) if s else np.nan) if t == 'N' else s)
    return {k: np.array(v) for k, v in out.items()}


def xyz_from_lonlat(lon, lat):
    """Unit vectors from longitude and latitude in degrees."""
    lo, la = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)
