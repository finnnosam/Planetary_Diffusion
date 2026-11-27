### 20 Nov 2026: a full-resolution surface map on demand (`scripts/render_surface.py`, `Render Surface.cmd`)

The user asked for an occasional full-resolution map of the surface at 0 Myr in the topography colours. The
snapshots are not reduced: every saved snapshot (`<seed>/snapshots/t+0000.0.npz`) holds every node at the plate
layer's resolution (n = 217, ~46 km per node at the equator); only the run's pictures are (`maps.png` samples a
360 × 180 grid at 70 dpi). The landscape stage's `present.npz` holds n = 868 (~11.5 km).

`scripts/render_surface.py` draws one surface at its grid's resolution: equirectangular, 4 n × oversample pixels
around the equator (default 2: 1,736 × 868 for the plate layer, 6,944 × 3,472 for the landscape stage), elevation in
`elevation_shading`'s fixed colours (65,536 levels, ±6,000 m) with its NW hillshade × 5 (the grid's spacing per row,
longitude wrapped), and a colour bar. Pixels between nodes are interpolated by inverse distance over the nearest node
and its eight neighbours (exact at nodes, never overshooting), so the shading has no nearest-node steps. It reads a
snapshot, a seed folder (`--time`, default 0), an arm or sweep folder (each arm's first seed) or a landscape run, and
writes beside the source (`<seed>/surface-t+0000.0.png`, `<run>/surface.png`); `--boundaries` overlays ridges,
trenches and collisions on a snapshot. Cloud timing: plate layer ~3 s; landscape stage 21 s, peak 1.9 GB (shaded in
row blocks, written with Pillow; a first version peaked at 5.3 GB). Checked on the seed-1020 pre-check snapshots and the
freeze-and-look run's 11.53 km surface (`outputs/planet-landscape/freeze-20261031-seed1001`).

`Render Surface.cmd` (root): drag a seed / arm / sweep folder, a snapshot or a landscape run onto it; double-click
renders the newest sweep's arms at 0 Myr. `tests/test_render_surface.py` (2; not in `PLANET_TESTS`: no model code).
Output only; nothing in the model or the runs changes.

Uncommitted on the device: `scripts/render_surface.py`, `tests/test_render_surface.py` (new). Suggested commit
subject: "Full-resolution surface render". Not published: `Render Surface.cmd`.
