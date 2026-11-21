"""Native-resolution map of a surface: elevation in the fixed topography colours (planet_reviews.elevation_shading,
65,536 levels) with a NW hillshade, on an equirectangular image at the resolution the model ran at.

The run's maps.png samples a 360 x 180 grid (~111 km) at 70 dpi; the snapshots hold every node of the grid the model
computes on (plate layer n = 217, ~46 km per node; the landscape stage's present.npz n = 868, ~11.5 km). This draws
one surface at that resolution, on demand:

  python -m scripts.render_surface outputs\\planet-calibration\\<sweep>\\<arm>\\seed1001          # its 0 Myr snapshot
  python -m scripts.render_surface outputs\\planet-calibration\\<sweep>                          # each arm's first seed
  python -m scripts.render_surface outputs\\planet-landscape\\<run>                              # the 11.53 km surface
  python -m scripts.render_surface <snapshot.npz> --time -100 --boundaries

A seed folder renders snapshots/t+0000.0.npz (or --time T). The image is written beside its source
(`<seed>/surface-t+0000.0.png`, `<run>/surface.png`): one pixel per node at the equator (the grid has 4 n nodes
there; 868 x 434 px for the plate layer, 3,472 x 1,736 for the landscape stage), plus a colour bar. Every pixel shows
its nearest node's own value: nothing is interpolated. A flat map cannot give every node exactly one pixel: toward the
poles neighbouring pixels repeat a node, and near the cube grid's face corners (mid-latitudes) the nodes are denser
than the pixels, so at native size ~6% of nodes have no pixel of their own (n = 217). --all-nodes samples a 1.5x finer
pixel grid instead, on which every node appears (still nearest-node values only, so some nodes cover two pixels).
--boundaries overlays the plate boundaries (ridges red, trenches black, collisions magenta) on a snapshot.
Elevations are relative to sea level.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere

from planet_tectonics.params import COLLISION, RIDGE, SUBDUCTION


def grid_for(size):
    """The CubeSphere whose node count is size."""
    n = int(round(np.sqrt(size/6)))
    for m in (n, n-1, n+1):
        g = CubeSphere(m)
        if g.size == size:
            return g
    raise SystemExit(f'no cube-sphere grid has {size} nodes')


def pixel_dirs(height, width, r0, r1):
    lat = np.pi/2-np.pi*(np.arange(r0, r1)+.5)/height
    lon = -np.pi+2*np.pi*(np.arange(width)+.5)/width
    la, lo = np.meshgrid(lat, lon, indexing='ij')
    return np.stack((np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)), -1).reshape(-1, 3)


def sample(g, z, height, width, rows=64):
    """Every pixel of the (height, width) equirectangular grid takes its nearest node: (values float32, node index)."""
    near = np.empty((height, width), np.int32)
    for r0 in range(0, height, rows):
        r1 = min(height, r0+rows)
        near[r0:r1] = fields.nearest_nodes(g, pixel_dirs(height, width, r0, r1)).reshape(r1-r0, width)
    return np.asarray(z, np.float32)[near], near


def shade(z, radius, exaggeration=5., rows=256):
    """elevation_shading.image's colours and NW light (uint8 RGB), with the equirectangular grid's spacing per row and
    the longitude wrapped; done in row blocks to keep memory near the image's own size."""
    from planet_reviews.elevation_shading import CMAP, NORM
    h, w = z.shape
    dy = radius*np.pi/h
    az, alt = np.radians([315, 45])
    light = np.array([np.sin(az)*np.cos(alt), np.cos(az)*np.cos(alt), np.sin(alt)], np.float32)
    out = np.empty((h, w, 3), np.uint8)
    for r0 in range(0, h, rows):
        r1 = min(h, r0+rows)
        a0, a1 = max(0, r0-1), min(h, r1+1)                           # one halo row each side
        zb = z[a0:a1]
        lat = np.pi/2-np.pi*(np.arange(a0, a1)+.5)/h
        dx = np.maximum(radius*np.cos(lat)*2*np.pi/w, dy*.05).astype(np.float32)[:, None]
        gx = (np.roll(zb, -1, 1)-np.roll(zb, 1, 1))/(2*dx)
        gy = np.gradient(zb, axis=0)/np.float32(dy)                   # rows run north to south
        k0, k1 = r0-a0, r0-a0+(r1-r0)
        nx, ny = -gx[k0:k1]*exaggeration, gy[k0:k1]*exaggeration
        inten = np.clip((nx*light[0]+ny*light[1]+light[2])/np.sqrt(nx**2+ny**2+1), 0, 1)
        rgb = CMAP(NORM(z[r0:r1]), bytes=False)[..., :3].astype(np.float32)
        out[r0:r1] = (np.clip(rgb*(.48+.65*inten[..., None]), 0, 1)*255+.5).astype(np.uint8)
    return out


def sources(path, time):
    """(npz file, output png, title) for a snapshot, a seed folder, an arm folder, a sweep folder or a landscape run."""
    p = Path(path)
    if p.is_file():
        return [(p, p.with_name(f'surface-{p.stem}.png'), str(p))]
    if (p/'present.npz').exists():
        return [(p/'present.npz', p/'surface.png', f'{p.name} (landscape stage)')]
    if (p/'snapshots').is_dir():
        f = p/'snapshots'/f't{time:+07.1f}.npz'
        if not f.exists():
            raise SystemExit(f'{f} not found (saved times: see {p/"snapshots"})')
        return [(f, p/f'surface-{f.stem}.png', f'{p.parent.name}/{p.name}, {time:g} Myr')]
    seeds = sorted(p.glob('seed*'))
    if seeds:
        return sources(seeds[0], time)
    arms = sorted(d for d in p.iterdir() if d.is_dir() and any(d.glob('seed*')))
    if arms:
        return [s for a in arms for s in sources(a, time)]
    raise SystemExit(f'{p}: not a snapshot, seed, arm, sweep or landscape folder')


def colour_bar(width, title, n, km, ppn, dpi=100, height=110):
    """The colour bar strip under the map, as uint8 RGB (width x height)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from planet_reviews.elevation_shading import CMAP, NORM
    fig = plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
    cax = fig.add_axes([.15, .62, .7, .14])
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation='horizontal')
    cb.set_ticks([-6000, -4000, -2000, 0, 1000, 2000, 4000, 6000])
    cb.ax.tick_params(labelsize=7)
    fig.text(.5, .2, f'{title}\nelevation relative to sea level (m; colours clipped at ±6,000 m); native grid n = {n}, '
             f'~{km:.0f} km per node at the equator, {ppn} px per node there; nearest-node values; NW hillshade ×5',
             ha='center', va='center', fontsize=6.5, linespacing=1.4)
    fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return img


def render(npz, png, title, all_nodes=False, boundaries=False):
    from PIL import Image
    d = np.load(npz)
    z = d['elevation_m'].astype(np.float64)
    g = grid_for(z.size)
    height = 2*g.n if not all_nodes else 3*g.n         # native: one pixel per node at the equator (4 n around it)
    width = 2*height
    zi, near = sample(g, z, height, width)
    rgb = shade(zi, g.radius)
    del zi
    if boundaries and 'boundary_type' in d.files:
        b = d['boundary_type'][near]
        for code, col in ((RIDGE, (255, 0, 0)), (SUBDUCTION, (0, 0, 0)), (COLLISION, (255, 0, 255))):
            rgb[b == code] = col
        del b
    del near
    bar = colour_bar(max(width, 868), title, g.n, 2*np.pi*g.radius/1e3/(4*g.n), '1.5' if all_nodes else '1')
    if bar.shape[1] > width:                           # a narrow map: centre it over the bar's width
        pad = np.full((height, bar.shape[1], 3), 255, np.uint8)
        o = (bar.shape[1]-width)//2
        pad[:, o:o+width] = rgb
        rgb, width = pad, bar.shape[1]
    Image.fromarray(np.concatenate([rgb, bar], 0)).save(png, optimize=False)
    return png, (width, height+bar.shape[0])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path', help='a snapshot .npz, a seed folder, an arm or sweep folder, or a landscape run folder')
    ap.add_argument('--time', type=float, default=0., help='snapshot time in Myr for seed folders (default 0)')
    ap.add_argument('--all-nodes', action='store_true', help='a 1.5x finer pixel grid on which every node appears '
                    '(default: native, one pixel per node at the equator)')
    ap.add_argument('--boundaries', action='store_true', help='overlay plate boundaries (snapshots only)')
    a = ap.parse_args(argv)
    for npz, png, title in sources(a.path, a.time):
        out, (w, h) = render(npz, png, title, a.all_nodes, a.boundaries)
        print(f'{out}  ({w} x {h} px)', flush=True)


if __name__ == '__main__':
    sys.exit(main())
