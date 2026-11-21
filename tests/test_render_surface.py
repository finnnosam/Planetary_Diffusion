"""scripts.render_surface (20 Nov, journal 127): a full-resolution elevation map of one surface, beside its source."""
import numpy as np
from PIL import Image

from planet_grid.sphere_graph import CubeSphere
from scripts import render_surface as R


def test_renders_a_snapshot_at_native_resolution(tmp_path):
    g = CubeSphere(240)                             # 4 n = 960 px: wider than the colour bar's minimum
    z = 3000*g.xyz[:, 2]-500                       # a smooth field: north high, south deep
    bt = np.zeros(g.size, np.int8); bt[::50] = R.RIDGE
    f = tmp_path/'t+0000.0.npz'
    np.savez(f, elevation_m=z.astype(np.float32), boundary_type=bt)
    (png, size), = [R.render(*s, boundaries=True) for s in R.sources(f, 0.)]
    im = np.asarray(Image.open(png))
    assert size[0] == 4*240 and im.shape[1] == 4*240 and im.shape[0] == size[1] == 2*240+110
    assert im[2].astype(int).mean() > im[2*240-3].astype(int).mean()   # high land light, deep sea dark
    (png2, size2), = [R.render(*s, all_nodes=True) for s in R.sources(f, 0.)]
    assert size2[0] == 6*240


def test_pixels_take_their_nearest_nodes_own_values():
    g = CubeSphere(16)
    z = np.random.default_rng(0).normal(size=g.size)
    v, near = R.sample(g, z, 32, 64)
    assert np.array_equal(v, z.astype(np.float32)[near])          # nothing interpolated
    assert near.min() >= 0 and near.max() < g.size
    assert R.grid_for(g.size).n == 16


def test_native_size_at_the_plate_grid_and_all_nodes_shows_every_node():
    g = CubeSphere(217)
    _, near = R.sample(g, np.zeros(g.size), 2*217, 4*217)
    assert len(np.unique(near[216])) == 868                        # the equator rows: one node per pixel
    _, near = R.sample(g, np.zeros(g.size), 3*217, 6*217)
    assert len(np.unique(near)) == g.size                          # --all-nodes: every node drawn
