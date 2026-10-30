# Contract: plate snapshots (S2 → S4, reviews)

Written by `scripts.simulate_plates` (one run) and `scripts.sweep_plates` (several parameter combinations × seeds),
from `planet_tectonics.PlateSim.state()` and `.plate_table()`. Read by `planet_landscape.coupling` (the landscape's
uplift and load forcing over the last 30–40 Myr) and by the reviews in `planet_reviews`.

## Folder layout

```
outputs/planet-plates/<YYYYmmdd-HHMMSS>/seed<s>/          one run (scripts.simulate_plates)
outputs/planet-calibration/sweep-<YYYYmmdd-HHMMSS>/       a sweep (scripts.sweep_plates)
    sweep.json             parameters shared by every combination
    results.json           the selection rule's result ('chosen', 'chosen_dir')
    <combination>/params.json  the parameters that differ (the folder is named by them)
    <combination>/seed<s>/ one run, laid out as below
    reviews/               review printouts

seed<s>/
    stats.json             time series and present-day statistics against Earth; 'params', 'params_effective'
    index.html, maps.png, present.png
    snapshots/t<time>.npz        world fields (below), time formatted as {t:+07.1f}: t-0250.0.npz … t+0000.0.npz
    snapshots/t<time>-plates.json  the plates at that time
```

Snapshots are saved every `--save-every` Myr (default 1) over the last `--save-last` Myr (default 40) and every
10 Myr before that. The landscape needs every 1 Myr over its window.

## `t<time>.npz`: world fields at one time

Node vectors on `CubeSphere(217)` (282,536 nodes), except the scalars.

| Field | dtype | Unit | Meaning |
|---|---|---|---|
| `time_myr` | f8 scalar | Myr | time of the snapshot (≤ 0) |
| `plate_id` | i4 | – | owning plate (ids match `-plates.json`) |
| `continental` | bool | – | continental crust (flag; the reviews also use a label-free test on thickness) |
| `crust_thickness_m` | f4 | m | crust thickness, sediment excluded |
| `crust_age_myr` | f4 | Myr | time since the crust formed |
| `ocean_age_myr` | f4 | Myr | ocean-floor age; NaN on continental crust |
| `orogeny_age_myr` | f4 | Myr | time since the last orogeny at the node |
| `boundary_type` | i1 | – | 0 none, 1 ridge, 2 subduction, 3 collision, 4 transform (`planet_tectonics.params`) |
| `convergence_mm_yr` | f4 | mm/yr | convergence rate across the boundary (0 off boundaries) |
| `elevation_m` | f4 | m | surface elevation relative to sea level (flexural isostasy, sediment, dynamic topography) |
| `sediment_m` | f4 | m | sediment: solid (grain) thickness with `sediment_compaction` (the default), else the layer's thickness |
| `sediment_bulk_m` | f4 | m | sediment column height with porosity (compaction) |
| `dynamic_topography_m` | f4 | m | mantle-driven surface motion, already added to `elevation_m` |
| `new_crust` | bool | – | nodes no plate covered in the last step, filled with new ocean crust |
| `sea_level_m` | f4 scalar | m | sea level on the fixed datum (fixed ocean-water volume) |
| `tectonic_thickening_m` | f4 | m | crust thickening by convergence and collision in the last step (before erosion) |

Runs made before 29 Oct 2026 also carry `extension_age_myr` and `thickening_age_myr` (from code since removed);
readers must ignore fields they do not know and must not require these two.

## `t<time>-plates.json`

A list, one entry per plate: `id`, `area_km2`, `continental_km2`, `speed_cm_yr` (at 90° from its Euler pole),
`axis` (unit Euler-pole vector, 3 floats) and `Q` (the plate's accumulated rotation, 3×3, row-major nested lists),
so a plate-fixed point p at that time sits at `Q @ p`.
