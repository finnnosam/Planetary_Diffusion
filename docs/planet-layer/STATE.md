# Planet layer: current state

**As of 27 Nov 2026, after the second independent review (journal 151; last pass 111, journal 150).** Rewrite this
file at the end of every pass so it is true; keep it under ~120 lines. The full STATE as it stood after pass 110 is in
`archive/STATE-2026-11-26-pass110.md`; the traps are in `gotchas.md`.

## Target

What the map shows against Earth on the same grid, relative to sea level: how much land, where it is (inland
lowlands, dry interiors, coastline roughness), how high it is (height distribution, mountains, plateaus). The
continents' level on the datum is a diagnostic, not the target: when the two disagree, the target decides and the
datum regression becomes the next problem (user, 26 Nov; journal 149; it replaces the 19 Nov target, journal 122 §6).
The 27 Nov review judged it sound; its wording was Claude's, and the user confirmed it on 27 Nov (journal 152).

## Right now

- **Nothing is running.** Focus stays on the plate stage (user, 27 Nov). The landscape and meso findings of the
  review are recorded for later (below; journal 151 §6).
- **Next (pass 112): repair C's return transfer.** C = `remap_rule` 'overlap' + `rift_area` 'plates' (built, off). It
  thins oceanic basement below 7 km with nothing booked, from two sources (`outputs/review-astra/2026-11-27/remap_trace.py`;
  seed 1001 over 250 Myr: 29,591 crossings):
  (a) local: forearc erosion cuts a world cell that combines several plate cells, and the return spreads the cut over
  every contributor, including cells with nothing to remove; this starts it (all of the first 10 Myr on 1001 and 1002;
  55% over the run);
  (b) the per-plate net correction (pass 109) moves untouched basement (45% over the run, 66% in the last 50 Myr).
  Order: a unit test first on a small physical example (a basement patch beside an accreted arc; a known forearc cut;
  debit only against available material; requested and removed volumes both conserved, any shortfall booked; zero
  change on distant basement; subaerial erosion apart; repeated with a newly exposed cell and an ownership change).
  **No global clamp** (it creates crust while passing the floor test). Then `remap_trace.py` (10 steps on 1001 and
  1002, 250 on 1001: no unbooked crossing), the planet tests with C as the default (pin journal 150 §3's three
  trajectory checks with notes), then re-sweep C against the control with a new frozen prediction
  (`precheck-p111-confirm/prediction.md` as template; it must also predict `shelf200` and `coast_step`), running
  2001–2006 in the same sweep.
- Next unused pre-check seed: 1033. The last default change to the plate layer: pass 104 (`rift_traction`).
- Independent review: done 27 Nov (`docs/review/2026-11-27-review.md`); next by about pass 126, or sooner if passes
  keep ending without progress on the target.

## The model

`PARAMS` in `packages/planet-tectonics/src/planet_tectonics/params.py` are the current model: a bare
`simulate_plates` run is it. n = 217 (46 km), 1 Myr steps, 250 Myr from a supercontinent. Where each value comes from:
`docs/planet-parameter-audit.md`. The settings by area and their dates: the archived STATE §3. Built but off: C
(conservative remap between the world grid and the plate rasters, the rift strained by the area the plates open);
not adopted in pass 111 (oceanic basement thinned unbooked).

## Scorecard against Earth (headline)

Six development seeds, mean over the saved snapshots −40 … 0 Myr (`geography`; the Rift Area Sweep's control
`sweep-20261126-161947/remap_rule-nearest_rift_area-world`, the model's six runs; the last four lines added
27 Nov, `reviews/geography-151.txt`), unless marked °: 0 Myr, six-seed medians from the 27 Nov review or pass 104's reviews.

| Measure | Model | Earth | |
|---|---|---|---|
| land, % of the globe | 27.3 | 28.1 | near |
| low flat interior (four dry rings, centre < 500 m, relief < 300 m), % of land | 3.8 | 14.2 | **main failure** |
| land below 200 m among land ≥ 300 km from water, % | 3.6 | 19.4 | **main failure** |
| same, distance from the main water body only, % | 11.9 | 22.0 | short |
| land ≥ 300 km from water / from the main water body, % of land | 39.1 / 53.0 | 62.1 / 68.6 | short |
| land–water crossings per Mkm²; one-cell holes / islands | 53.4; 358 / 350 | 19.1; 13 / 29 | rough |
| shore step (median height jump across a land–water edge) | 303 m | 132 m | steep |
| water 0–200 m deep, % of water | 6.4 | 6.8 | near |
| land p50 / p99, m ° | 444 / 3,066 | 379 / 4,680 | middle high, top short |
| land above 4 km, % of land; largest connected area above 4 km ° | 0.15; 0.08 Mkm² | 1.90; 2.10 Mkm² | no plateaus |
| median height of continent at 60–65 km total crust, 1° cells ° | 2,855 m | 4,023 m | thick crust stands low |
| crust ≥ 30 km below the datum, share (diagnostic) ° | 0.193 | 0.112 | low tail |
| continental crust area / median thickness ° | 39.5% / 37.3 km | 38.4% / 37.0 km | good |
| ocean floor production, last 100 Myr ° ; floor < 10 Myr ° | 5.19 km²/yr; 47.4 Mkm² | 3.37; 29.9 | renews too fast |
| deep-floor neighbour difference ° | 39 m | 77 m | too smooth |
| crystalline continent lost, elapsed 250 → 500 Myr (seeds 1001 / 1002) ° | 16.7 / 10.7% | — | late drift |

C on the same seeds (window; not adopted): holes 358 → 213, crossings 53.4 → 43.0, low flat interior 3.8 → 6.8%, far
low land 3.6 → 9.1% (main water body 11.9 → 16.6%), far land 39 → 46%, land 27.3 → 27.5%; its costs: water 0–200 m
6.4 → 3.8% (6 of 6 seeds), shore step 303 → 365 m, continents lower on the datum (0.193 → 0.281).

## Open problems (plate stage), in order

1. **Rough coastlines and missing inland lowlands.** Traced to the composition step (journals 143–147); the candidate is
   C, blocked by its unbooked thinning of basement (pass 112, above). Read inland distance with both water definitions
   side by side: the difference does not say which depressions should be lakes, dry basins or seas.
2. **The lithosphere under C, in km and m** (journal 151 §2, replacing the "thermal-age trace"). C's continents sit
   lower on the datum. Effective age (5.5 Myr per Myr under C) is an inverse-table clock and cannot measure it. On
   persistent parcels, split the change in L, column mass and uplift into mixing, mechanical strain, cooling and
   resets, both remaps, with the rule written first.
3. **Plateaus: thick crust stands low and the thickest roots are missing** (unparked 27 Nov). There is more 60–65 km
   crust than on Earth, 1.2 km lower; there are no 70–80 km cells (Earth 0.81 Mkm²); 45–50 km crust is already higher
   than Earth's, so no datum rise and no general thickening. Pick the broadest sustained collision by its incoming
   volume and duration before reading elevation. Book flux, root growth, spreading, foundering (65 km, τ 10 Myr) and
   the 'pressure' density (an upper bound). Not a Tibet per seed.
4. **The ocean basin.** Too much young floor while ridge length matches Earth's. First integrate normal opening × ridge
   length against the created area and the age cohorts. Deep floor too smooth within age: separate, and no further
   age–depth fit. All ocean–ocean convergence is trench (instant initiation).
5. **Late crust loss in long runs**: a lifetime constraint of its own. Book collision, rifting and sediment-to-crust
   conversions by interval on the two 500 Myr series before any longer run; then a supercontinent cycle.
6. Parked: the dynamic topography's placement; the passive edge (0.56× Earth's length) and the breakup lottery.

## Later stages (recorded, not worked on now; journal 151 §6, `planet-ideal-flow.md` S4 / S7)

- **Landscape: the largest loss in the finished terrain.** On seed 1001's chain, places above 3 km on the plate
  parent stand at 4,448 m with forcing alone and at 638 m finished. On Earth the law overpredicts basin denudation
  11.4× (8.4× without the slope cap). The present erodibility is applied through all 30 Myr. First, when the stage
  resumes: the incision law on Earth (the spatial support of the 11.5 km cell; residuals by runoff and steepness),
  then a replay of the belts at parcel-time erodibility on two seeds.
- **Meso:** to be redone; freeze Earth's test regions, classes and preprocessing first. Earth-input and planet ratios
  do not split the gap between causes (journal 146 corrected).

Corrections queued from the 22 Nov review: `gotchas.md`, last section.

## Where things are

| What | Where |
|---|---|
| working rules | `CLAUDE.md` (repo root) |
| decisions (one row each), the journal (one file per pass), traps | `decisions.md` and archives, `journal/`, `gotchas.md` |
| code, scripts, data, commands | `docs/project-map.md` |
| the independent reviews (22 and 27 Nov) and their probes | `docs/review/`, `outputs/review-astra/` (`2026-11-27/` for the second) |
| probes of each pass | `outputs/planet-calibration/precheck-pNNN-*/` (the pass's journal names its folder) |
| sweeps and their reviews | `outputs/planet-calibration/sweep-*/`, `reviews/` inside; one root `.cmd` each |
