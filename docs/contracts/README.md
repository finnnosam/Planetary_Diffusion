# Stage contracts

What each stage writes and what the next stage reads. A stage can be rewritten freely as long as its outputs keep
to its contract; a change to a contract is a change to every stage that reads it, so it is made here first.

| Contract | Written by | Read by | Status |
|---|---|---|---|
| [Plate snapshots](plate-snapshots.md) | `planet_tectonics` (S2, via `scripts.simulate_plates` / `scripts.sweep_plates`) | `planet_landscape.coupling` (S4), `planet_reviews` | in use |
| [23 km hand-off](handoff-23km.md) | `planet_landscape` (S4, via `scripts.run_planet_landscape`) | `planet_meso` (S7), `scripts.planet_meso_check` | in use (elevation only) |
| [Climate](climate.md) | `planet_climate` (S3, S5) | S4, S6 | draft: nothing is built yet |

Stage numbers (S0–S9) are those of [`planet-ideal-flow.md`](../planet-ideal-flow.md).

## Shared conventions

- **Grid.** Every planet field lives on the shared-node equi-angular cube of `planet_grid` (`CubeSphere(n)`), as a
  node vector of length 6 n² + 2 in `CubeSphere` node order. `g.to_cube(v)` gives the (6, n+1, n+1) face charts;
  `g.xyz` the unit vectors; `g.area` the node cell areas (m², they tile the sphere). Radius 6,371 km unless stated.
  Grids in use: n = 217 (46 km, plates), n = 434 (23 km, hand-off), n = 868 (11.53 km, landscape),
  N = 3,472 (2.88 km, meso).
- **Units.** Metres, years (time in Myr, negative before the present: −250 … 0), m²/yr, mm/yr; float32 fields,
  int32 ids, bool masks.
- **Elevation** is relative to the run's sea level (0 m = shoreline); the sea level on the fixed datum is stored
  beside it where it moves.
- **Files.** `.npz` (numpy, compressed) for fields, `.json` for tables and parameters. Every run saves the
  parameters it actually used (`params_effective`), so older runs stay readable when defaults change.
