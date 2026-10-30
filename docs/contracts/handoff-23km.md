# Contract: the 23 km hand-off (S4 → S7)

Written by `scripts.run_planet_landscape` (the landscape stage, `planet_landscape`). Read by the meso stage
(`planet_meso`) and checked by `scripts.planet_meso_check`. Today the hand-off carries elevation only; the full
conditioning bundle (S6 in [`planet-ideal-flow.md`](../planet-ideal-flow.md): land fraction, rivers, lakes, roughness
type, grain, lithology, climate) adds channels to this contract when each is built.

## Landscape output: `present.npz`

`outputs/planet-landscape/<YYYYmmdd-HHMMSS>-seed<s>/present.npz` (or `theta-*/theta<k>/present.npz` from the theta
calibration), with `stats.json` and `preview.png` beside it. Node vectors on `CubeSphere(n)`, n = 868 (11.53 km,
4,520,546 nodes):

| Field | dtype | Unit | Meaning |
|---|---|---|---|
| `n` | int scalar | – | cube intervals per face (868) |
| `elevation_m` | f4 | m | present-day surface relative to sea level |
| `drainage_area_m2` | f4 | m² | upstream drainage area |
| `outlet` | i8 | – | node each node finally drains to |
| `river` | bool | – | land with drainage area ≥ 2·10¹⁰ m² |
| `lake_depth_m` | f4 | m | depth of the filled surface above the ground (lakes, closed basins) |
| `ocean` | bool | – | below sea level and connected to the ocean |
| `drains_to_sea` | bool | – | land whose outlet is ocean |
| `rock_load_m` | f4 | m | net rock added (uplift minus erosion plus deposition) over the window |
| `erodibility` | f4 | – | erodibility factor used |
| `initial_m` | f4 | m | the surface at the start of the window |
| `denudation_late_m_myr` | f4 | m/Myr | net rock removal over the last 10 Myr |

## The 23 km input of the meso model

The meso model reads 23.04 km **cells** of its own cube charts (n = 434, six faces of 434 × 434 cells), not nodes.
Because the landscape grid is exactly twice as fine (868 = 2 × 434), the hand-off needs no resampling between grids:

1. `fine = CubeSphere(868).to_cube(present['elevation_m'])`, the (6, 869, 869) face charts;
2. each 23 km cell (face f, row r, column c) is the mean of 8 × 8 points of `fine`, read bilinearly with
   `planet_grid.cube.sample` at the chart positions of the fine nodes 8c … 8c+7 (cell centres at c + 7/16),
   exactly as `planet_meso.planet_tiled` prepares its own input (`scripts.planet_meso_check.cells23`);
3. the model sees the cells as signed-square-root metres, `sign(z) · sqrt(|z|)` (`planet_meso.sphere.encode`),
   in windows of 18 × 18 cells (16 × 16 core, 368.64 km) cut from the face charts.

The meso's output is 2.88 km elevation on the N = 3,472 cube (8 × finer than the input cells), in metres.

**Rule for every future channel** (cross-cutting rule 1 of the ideal flow): a channel enters this contract only with
the same operator applied on Earth for training and on the planet for inference.
