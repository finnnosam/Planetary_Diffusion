# Speeding up plate tectonics sweeps

Start with the cube-coordinate lookup and the mantle velocity calculation. Both run on the same 282,536-node grid throughout a 250 Myr trajectory. The first can answer most coordinate queries with a small integer table. The second can reuse the grid geometry and skip a full velocity calculation on the default path. Implement and time each independently so its contribution is visible.

The measurements below are isolated benchmarks from 27 November 2026. Unless specified otherwise, they used seed 1001 snapshots from the adopted overlap-remap arm at production resolution. Sweep wall time is the deciding measure. The previous speed pass reduced single-run step time by 14% at two threads and concurrent-run step time by 9%. [Previous profile and timings](journal/154_2026-11-27_speed.md)

## First implementation pass

### 1. Use a small table for cube coordinates

For each rotated point, the current kernel picks a cube face and calls `atan2` twice to obtain two rounded chart indices. Keep the face selection and tie rules. On the selected face, both coordinate ratios lie in the interval from −1 to 1. Divide that interval into bins. Store the integer chart index when a conservatively widened bin lies wholly within one rounded cell. Mark bins touching a rounding boundary so that queries in those bins execute the existing formula.

Build the table once per grid resolution and reuse it for both axes and every plate. The prototype used 16,384 bins. Across 66 plate rotations at three epochs, all 18.65 million chart entries matched. Median time per full-grid rotation fell from 2.73–2.81 ms to 0.93–0.96 ms on two threads. Before enabling it, exercise both sides of half-cell thresholds, cube edges and corners, face-selection ties and the nine-point overlap queries. [Lookup kernel](../../packages/planet-grid/src/planet_grid/fields.py) · [Overlap remap](../../packages/planet-tectonics/src/planet_tectonics/composition.py)

### 2. Reuse mantle grid geometry and calculate velocity on demand

Compute the harmonic transform coordinates and the two local tangent vectors once on the fixed world grid. Pass them into scalar synthesis and gradient evaluation. Keep the coordinate calculation for moving slab parcels and arbitrary evaluation points. On one saved velocity field, caching the geometry reduced median gradient time from 30.9 to 16.0 ms with identical output, at a cost of 18.1 MB per worker.

Each flow update also constructs an **unscaled** velocity field. The default has basal traction off. Rift traction reads a separately scaled field, and snapshots save velocity coefficients. Materialize the unscaled field when a caller requests it, invalidate it whenever the coefficients change, and retain the existing initialization and direct-access behavior. Exercise basal traction enabled, measured-power rift traction and the unscaled flow mode. Skipping one evaluation saves approximately the cost of one gradient on steps where it goes unread. Count that saving using the version of the gradient actually in use. [Mantle update and gradients](../../packages/planet-tectonics/src/planet_tectonics/mantle_flow.py) · [Motion reader](../../packages/planet-tectonics/src/planet_tectonics/motion.py) · [Rift traction tests](../../tests/test_rift_traction.py)

## Additional measured opportunities

| Change | Implementation and evidence |
|---|---|
| Fill gaps from their borders | Specialize the ownership search in `_world_from`. For each unclaimed node, initialize from adjacent claimed sources using the original directed edge length and source-index tie rule. Propagate into the remaining gaps with the existing distance and node ordering. Across four saved gap masks, only 2,283–3,483 cells needed labels. Labels and distances matched, and the call took 0.86–1.16 ms versus 27–29 ms. Include empty gaps, multi-cell gaps and equal-distance ties in the comparison. [Ownership](../../packages/planet-tectonics/src/planet_tectonics/composition.py) · [Reference search](../../packages/planet-tectonics/src/planet_tectonics/graphops.py) |
| Precompute shelf edges | Store the 565,068 undirected axis edges in their current traversal order, with endpoints, distance and width. Rebuild the active mask from the current ocean and depth on every shelf call. On three saved epochs, the 13.6 MB edge table saved 5–8 ms per call and matched elevation, available material and volume moved off land. Four calls per Myr imply roughly 5–8 seconds per 250 Myr run at those rates. [Shelf transport](../../packages/planet-tectonics/src/planet_tectonics/shelf.py) |
| Persist the crust-flow tables | Sheet and channel rheology are already tabulated in each process. Store those exact arrays for reuse across workers and sweeps. Include physical parameters, table axes, integration resolution and implementation version in the key. Publish a complete file atomically so concurrent workers can read it safely. A fresh two-thread process took 0.24 s for the sheet table and 6.35 s for the channel table. The pair occupies 1.82 MiB. This saves worker startup time. [Table builder and current cache](../../packages/planet-tectonics/src/planet_tectonics/crust_flow.py) |

## Candidates to benchmark next

- **Reuse return-transfer geometry for crust, sediment and rock.** Within one plate update, compute contributor weights, destination weight totals and reached-node indices once. Keep each material’s availability, floor, gain and loss accounting, and the current update order. The earlier profile attributed 6% of a step to the whole return operation. Only the shared portion is recoverable. [Material return](../../packages/planet-tectonics/src/planet_tectonics/composition.py)
- **Cache one interpolation coordinate during crust flow.** Lithosphere thickness stays fixed throughout one call to `flow`. Precompute its edge interpolation index and weight, then reuse them across substeps as crust thickness and slope change. Reusing the *entire* previous edge result had little promise in three early production-grid steps on seed 9999. The mean fraction of edges with both thickness inputs unchanged was about 0.2% on step one and below 0.1% on each of the next two. [Crust-flow loop](../../packages/planet-tectonics/src/planet_tectonics/crust_flow.py)
- **Reuse work in the per-Myr recorder.** Pass the sediment bulk thickness already calculated for elevation through to the same state record. Compute the plate table once for the summary and any snapshot at that epoch. Keep the existing timing of elevation evaluation because it also updates sea level. [State output](../../packages/planet-tectonics/src/planet_tectonics/sim.py) · [Elevation](../../packages/planet-tectonics/src/planet_tectonics/surface.py) · [Recorder](../../scripts/simulate_plates.py)

## Reduce elapsed time around the simulation

Use the existing borrowed-control option on confirmation sweeps when a matching confirmation reference exists. Retain the check-seed comparison and fallback. The Transfer development sweep ran seven trajectories in 24:05 after borrowing five controls. Its confirmation ran twelve in 43:31 on different seeds. Seven jobs fit in one wave on eight workers. [Transfer launcher](../../Transfer%20Sweep.cmd) · [Borrowing logic](../../scripts/sweep_plates.py)

After the simulations finish, try two or three concurrent workers for independent reviews. Pin each review to the intended sweep directory, check file dependencies, and retain individual exit status and logs. Reviews took 3:56 on that development sweep. Measure elapsed time and disk contention before increasing concurrency. [Review runner](../../scripts/run_sweep.py)

## Verify identity and measure the sweep

1. Capture a baseline from the current checkout. Repeat the repository’s fingerprint cases on each change independently: grid 48, seed 7, 70 steps at one and two threads, then grid 217, seed 9999, 60 steps at two threads. Require every recorded step to match. For changes to lazy fields or caches, also compare the returned values on first access and after invalidation. [Fingerprint script](../../outputs/planet-calibration/precheck-speed/ident.py) · [Comparison script](../../outputs/planet-calibration/precheck-speed/compare.py)
2. Time the same workload before and after. Separate cold table construction from warm per-step time. Include a late epoch with many plates. Then measure eight concurrent runs at two threads each and record peak memory. The shelf and mantle caches together add about 32 MB per worker before any sharing. [Timing harness](../../outputs/planet-calibration/precheck-speed/timing.py)
3. Run the usual sweep tests and control identity check. Record simulation, tests, reviews and total elapsed time separately. Use the observed end-to-end result when deciding which changes to keep. [Sweep workflow](../../scripts/run_sweep.py)

## Where a learned model could help

A learned model is a separate experiment for screening candidate settings. First test how well a cheap predictor ranks changes in land fraction, inland lowlands and coastline roughness from development runs. A conditional diffusion model could then be tested on the same task using coarse plate states to generate finer elevation detail. Compare ranking accuracy, training cost and inference time on entire held-out development trajectories and settings. Keep seeds 2001–2006 reserved for physical confirmation. This experiment becomes worthwhile if it reduces the number of full runs needed to identify promising candidates. [Project targets and seed rules](../../CLAUDE.md) · [Fields available for training](../contracts/plate-snapshots.md)
