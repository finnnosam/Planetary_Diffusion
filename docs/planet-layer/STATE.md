# Planet layer: current state

**As of 27 Nov 2026, after pass 112 (journal 153).** Rewrite this file at the end of every pass so it is true; keep it
under ~120 lines. The full STATE as it stood after pass 110 is in `archive/STATE-2026-11-26-pass110.md`; the traps are
in `gotchas.md`.

## Target

What the map shows against Earth on the same grid, relative to sea level: how much land, where it is (inland
lowlands, dry interiors, coastline roughness), how high it is (height distribution, mountains, plateaus). The
continents' level on the datum is a diagnostic, not the target: when the two disagree, the target decides and the
datum regression becomes the next problem (user, 26 Nov; journal 149; confirmed by the user 27 Nov, journal 152).

## Right now

- **Nothing is running.** Focus stays on the plate stage (user, 27 Nov).
- **Pass 112 adopted C** (`remap_rule` 'overlap' + `rift_area` 'plates', with the repaired material-aware return) as a
  numerics fix: every veto on 1001–1006 and 2001–2006, oceanic basement never below 7 km (Transfer Sweep:
  `sweep-20261127-192554` development, `-195437` confirmation). **Its frozen terrain prediction failed on the land
  fraction** (−2.39 points on 2001–2006, band −1.5 … +2.0; development −0.87), so no terrain claim. The other six
  deciding lines held.
- **Next (pass 113): why the continents stand lower under C, in km and m**, aimed at the land loss. On 2001–2006 the
  share of crust ≥ 30 km under the datum goes 0.25 → 0.34 (Earth 0.11), and `three_factor`'s level term on the land
  median goes −12 → −68 m. On persistent continental parcels, split the change in lithosphere thickness L, column mass
  and uplift into mixing, mechanical strain, cooling and resets, under C and under 'nearest' / 'world' (now the
  comparison arm), with the rule written first. Effective age is reported only as a symptom (it is an inverse-table
  clock, journal 151 §2). The candidate fixes (thermal age read at the nearest plate node under 'overlap'; thickening
  that keeps the age on Lamb's flat stretches; journal 150 §5) wait for that split.
- Next unused pre-check seed: 1033 (pass 112 used none). Last default change to the plate layer: pass 112 (C).
- Independent review: done 27 Nov (`docs/review/2026-11-27-review.md`); next by about pass 126, or sooner if passes
  keep ending without progress on the target.

## The model

`PARAMS` in `packages/planet-tectonics/src/planet_tectonics/params.py` are the current model: a bare
`simulate_plates` run is it. n = 217 (46 km), 1 Myr steps, 250 Myr from a supercontinent. Where each value comes from:
`docs/planet-parameter-audit.md`. The settings by area and their dates: the archived STATE §3. Since pass 112 the
world grid reads the plate rasters by first-order conservative remapping ('overlap'). The return is material-aware:
a loss is taken only from material above each node's floor, gains and losses are closed per plate without a uniform
offset, and any shortfall is booked. The rift strains by the area the plates open. Old runs read with
`crust_review.PREDATING` 'nearest' / 'world'.

## Scorecard against Earth (headline)

Window rows: six development seeds, mean over the saved snapshots −40 … 0 Myr, **the new default C**
(`sweep-20261127-192554/remap_rule-overlap_rift_area-plates`, `reviews/geography.txt`); in brackets the old default
A on the same seeds. Rows marked ° are 0 Myr six-seed medians from the 27 Nov review **on A**, not yet re-read on C.

| Measure | Model | Earth | |
|---|---|---|---|
| land, % of the globe | 26.4 (27.3) | 28.1 | **short, worse under C** (2001–2006: 25.5 against 27.9) |
| low flat interior (four dry rings, centre < 500 m, relief < 300 m), % of land | 6.9 (3.8) | 14.2 | **main failure** |
| land below 200 m among land ≥ 300 km from water, % | 9.9 (3.6) | 19.4 | **main failure** |
| same, distance from the main water body only, % | 16.8 (11.9) | 22.0 | short |
| land ≥ 300 km from water / from the main water body, % of land | 45.5 / 56.4 (39.1 / 53.0) | 62.1 / 68.6 | short |
| land–water crossings per Mkm²; one-cell holes / islands | 42.6; 186 / 272 (53.4; 358 / 350) | 19.1; 13 / 29 | rough |
| shore step (median height jump across a land–water edge) | 368 m (303) | 132 m | steep |
| water 0–200 m deep, % of water | 3.7 (6.4) | 6.8 | **short under C** |
| land p50 / p99, m ° | 444 / 3,066 | 379 / 4,680 | middle high, top short |
| land above 4 km, % of land; largest connected area above 4 km ° | 0.15; 0.08 Mkm² | 1.90; 2.10 Mkm² | no plateaus |
| median height of continent at 60–65 km total crust, 1° cells ° | 2,855 m | 4,023 m | thick crust stands low |
| crust ≥ 30 km below the datum, share (diagnostic; C at 0 Myr) | 0.291 (0.193) | 0.112 | low tail, worse under C |
| continental crust area / median thickness ° | 39.5% / 37.3 km | 38.4% / 37.0 km | good |
| ocean floor production, last 100 Myr (C) ; floor < 10 Myr ° | 5.30 (5.19) km²/yr; 47.4 Mkm² | 3.37; 29.9 | renews too fast |
| deep-floor neighbour difference ° | 39 m | 77 m | too smooth |
| crystalline continent lost, elapsed 250 → 500 Myr (seeds 1001 / 1002) ° | 16.7 / 10.7% | — | late drift |

## Open problems (plate stage), in order

1. **Less land under C, and the continents lower on the datum** (pass 113, above). The land fraction is a target
   measure, and C took it away from Earth on the confirmation seeds. The cause is suspected in the lithosphere
   (journal 150 §5: C's L thicker in every crust bin), not yet measured in km and m.
2. **Rough coastlines and missing inland lowlands.** C moved every "where it is" measure toward Earth (holes −49%,
   crossings −19%, low flat interior +3.8 points, far low land +6.6 on 2001–2006), but the gaps stay large (above).
   C's shoreline costs: shallow water 6.4 → 3.7% of water, shore step +64 m. Read inland distance with both water
   definitions side by side.
3. **Plateaus: thick crust stands low and the thickest roots are missing.** There is more 60–65 km crust than on Earth,
   standing 1.2 km lower; there are no 70–80 km cells (Earth 0.81 Mkm²); 45–50 km crust already stands higher than
   Earth's. Pick the broadest sustained collision by its incoming volume and duration before reading elevation. Book
   flux, root growth, spreading, foundering (65 km, τ 10 Myr) and the 'pressure' density (an upper bound). It shares
   the vertical-support question with item 1. Not a Tibet per seed.
4. **The ocean basin.** Too much young floor while ridge length matches Earth's. Integrate normal opening × ridge
   length against the created area and the age cohorts first. The deep floor is too smooth within age: a separate
   problem, with no further age–depth fit. All ocean–ocean convergence is trench (instant initiation).
5. **Late crust loss in long runs**: a lifetime constraint of its own. Book collision, rifting and sediment-to-crust
   conversions by interval on the two 500 Myr series (both run on the old default A) before any longer run.
6. Parked: the dynamic topography's placement; the passive edge (0.56× Earth's length) and the breakup lottery; the
   plate pieces' T2 / T3 / T4 marked FAIL under C (scored).

## Later stages (recorded, not worked on now; journal 151 §6, `planet-ideal-flow.md` S4 / S7)

- **Landscape: the largest loss in the finished terrain.** On seed 1001's chain, places above 3 km on the plate
  parent stand at 4,448 m with forcing alone and 638 m finished. On Earth the law overpredicts basin denudation 11.4×
  (8.4× without the slope cap), and the present erodibility is applied through all 30 Myr. First, when the stage
  resumes: test the incision law on Earth (the 11.5 km cell's spatial support; residuals by runoff and steepness), then
  replay the belts at parcel-time erodibility on two seeds.
- **Meso:** to be redone; freeze Earth's test regions, classes and preprocessing first.

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
