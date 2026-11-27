### 20 Nov 2026: the surface map at the model's native resolution (corrects journal 127)

The user wants the map at the resolution the model runs at, not enlarged. Journal 127's version drew 2 pixels per
node (`--oversample`, default 2) and filled pixels between nodes by inverse-distance blending: a larger, smoother
picture with values the model never computed. Changed in `scripts/render_surface.py`:

- **Native by default:** one pixel per node at the equator (4 n × 2 n: 868 × 434 map pixels for the plate layer,
  3,472 × 1,736 for the landscape stage), every pixel its nearest node's own value; nothing is interpolated.
  `--oversample` and the blending are removed.
- A flat map cannot give every node exactly one pixel. Measured at n = 217: the equator rows use 868 distinct nodes
  for 868 pixels, but the whole native image draws 93.8% of the nodes (the cube grid is denser than the pixels near
  its face corners, at mid-latitudes; toward the poles pixels repeat nodes). `--all-nodes` samples a 1.5× pixel grid
  (1,302 × 651), on which every node appears (1.25×: 99.6%), still nearest-node values only.
- The colour bar's label is two lines (it overflowed at 868 px). Landscape stage: 9 s in the cloud.
- `tests/test_render_surface.py` (3): native size and `--all-nodes` size; pixels equal their nearest node's value;
  at n = 217 the equator rows hold 868 distinct nodes and `--all-nodes` draws every node.

Uncommitted on the device: `scripts/render_surface.py`, `tests/test_render_surface.py` (new). Suggested commit
subject: "Native-resolution surface render".
