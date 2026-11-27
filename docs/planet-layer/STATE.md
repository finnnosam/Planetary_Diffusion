# Planet layer: current state

**As of 26 Nov 2026, after pass 111 (journal 150).** Rewrite this file at the end of every
pass so it is true; keep it under ~120 lines. The full STATE as it stood after pass 110 (model table, the long
scorecard, every open item) is in `archive/STATE-2026-11-26-pass110.md`; the traps it listed are in `gotchas.md`.

## Target

What the map shows against Earth on the same grid, relative to sea level: how much land, where it is (inland
lowlands, dry interiors, coastline roughness), how high it is (height distribution, mountains, plateaus). The
continents' level on the datum is a diagnostic, not the target: when the two disagree, the target decides and the
datum regression becomes the next problem (user, 26 Nov; journal 149; it replaces the 19 Nov target, journal 122 §6).

## Right now

- **Nothing is running.** Pass 111 read the Rift Area Sweep (`sweep-20261126-161947`): C (`remap_rule` 'overlap' +
  `rift_area` 'plates') passed every veto as the sweep ran them, but with C made the default the planet tests fail:
  the overlap remap thins oceanic crust below 7 km with nothing booked (1.9% of the ocean, 0.26 Mkm³ by 0 Myr at
  n = 217; also in the Remap Sweep's overlap arm, so it is the remap's). **C is not adopted**; both settings stay off.
- **Next (pass 112): the failing veto.** Trace where the overlap return thins oceanic crust (rule first; suspect: the
  return's per-plate correction of each crust type's net volume, pass 109, spread over all of a plate's oceanic
  nodes). Fix it; run the planet tests with C as the default (they must pass; pin the three trajectory checks of
  journal 150 §3 with notes); then re-sweep C against the control on 1001–1006, with a new terrain prediction frozen
  first so 2001–2006 run in the same sweep (`precheck-p111-confirm/prediction.md` is the template; rewrite its ranges).
- **Then: the lithosphere age.** Under C the continents' mean lithosphere age rises 5.5 Myr per Myr (median 3,600 Myr
  by 0 Myr; 1.3 per Myr under A). That is why they sit lower on the datum (share under it 0.193 → 0.281). Causes:
  Lamb's flat stretches turn any collision thickening into an age jump (both arms), and the overlap read mixes ages.
  Next arm: a per-step trace splitting the age change into aging, resets, thickening jumps and write-back (rule first).
- Next unused pre-check seed: 1033. The last default change to the plate layer: pass 104 (`rift_traction`).
- Due: an independent review (fresh subagent, code and data before docs) by about pass 125, or sooner: passes 105–111
  changed no default, each blocked by a new integrity problem while the geography under the remap keeps improving.

## The model

`PARAMS` in `packages/planet-tectonics/src/planet_tectonics/params.py` are the current model: a bare
`simulate_plates` run is it. n = 217 (46 km), 1 Myr steps, 250 Myr from a supercontinent. Where each value comes from:
`docs/planet-parameter-audit.md`. The settings by area and their dates: the archived STATE §3. Built but off: `remap_rule`
'overlap' (conservative remap between the world grid and the plate rasters) with `rift_area` 'plates' (the rift strains
by the area the plates actually open); not adopted in pass 111 (oceanic crust thinned unbooked).

## Scorecard against Earth (headline)

Today's model, six development seeds, mean over the saved snapshots −40 … 0 Myr (`geography`; the Remap Sweep's
control, `sweep-20261126-141105/remap_rule-nearest`), unless marked °: 0 Myr, `sweep-20261122-150810/combo-01`, from
the independent review or pass 104's reviews.

| Measure | Model | Earth | |
|---|---|---|---|
| land, % of the globe | 27.3 | 28.1 | near |
| low flat interior (four dry rings, centre < 500 m, relief < 300 m), % of land | 3.8 | 14.2 | **main failure** |
| land below 200 m among land ≥ 300 km from water, % | 3.6 | 19.4 | **main failure** |
| land ≥ 300 km from water, % of land | 39.1 | 62.1 | short |
| land–water crossings per Mkm² | 53.4 | 19.1 | 3× too rough |
| one-cell water holes / islands per globe | 358 / 350 | 13 / 29 | rough |
| land median / 99th percentile, m ° | 480 / 3,090 | 379 / 4,680 | middle high, top short |
| land above 4 km, % of land; largest connected area above 4 km ° | 0.19; 0.08 Mkm² | 1.90; 2.10 Mkm² | no plateaus |
| crust ≥ 30 km below the datum, share (diagnostic) ° | 0.193 | 0.112 | low tail |
| continental crust area / median thickness ° | 39.5% / 37.3 km | 38.4% / 37.0 km | good |
| ocean floor production, last 100 Myr ° | 5.19 km²/yr | 3.37 | renews too fast |
| relief within age on deep floor, 20–60 Myr ° | 127 m | 281 m | too smooth |
| crystalline continent lost, elapsed 250 → 500 Myr (seeds 1001 / 1002) ° | 16.7 / 10.7% | — | late drift |

On the development seeds C (overlap + plates, Rift Area Sweep, not adopted) moved holes 358 → 213, crossings
53.4 → 43.0, low flat interior 3.8 → 6.8%, far low land 3.6 → 9.1%, far land 39 → 46%, land 27.3 → 27.5%.

## Open problems, in order

1. **Rough coastlines and missing inland lowlands.** Traced to the composition step: each step the plate rasters,
   read back onto the world grid by nearest node, add ~160 one-cell holes; erosion repairs most of it (journals
   143–147). The candidate is the conservative remap (above), blocked by its thinning of oceanic crust (pass 111).
   After it: the height distribution of the land.
2. **The landscape stage's missing relief.** The orogeny-age hand-off was fixed in pass 105; the landscape has no
   channel for the lithosphere's support and its erosion law's form is not Earth's (journals 142, 144). Behind item 1.
3. **The ocean basin.** Too much young floor while ridge length matches Earth's; all ocean–ocean convergence is trench
   (instant initiation); deep floor too smooth within age. Integrate opening along ridges before touching their length.
4. **Late crust loss in long runs** (collision consumes area faster than rifting restores it). Trace the crust budget's
   terms every 50 Myr on the 500 Myr runs. Then a supercontinent cycle for longer runs (`docs/planet-ideal-flow.md`).
5. Parked: the dynamic topography's placement; broad high plateaus (test the capacity given a suitable collision, not a
   Tibet per seed); the passive edge (0.56× Earth's length) and the breakup lottery; the meso stage (to be redone with
   more inputs after the geography settles).

Corrections queued from the independent review: `gotchas.md`, last section.

## Where things are

| What | Where |
|---|---|
| working rules | `CLAUDE.md` (repo root) |
| decisions (one row each), the journal (one file per pass), traps | `decisions.md` and archives, `journal/`, `gotchas.md` |
| code, scripts, data, commands | `docs/project-map.md` |
| the independent review (22 Nov) and its probes | `docs/review/2026-11-review.md`, `outputs/review-astra/` |
| probes of each pass | `outputs/planet-calibration/precheck-pNNN-*/` (the pass's journal names its folder) |
| sweeps and their reviews | `outputs/planet-calibration/sweep-*/`, `reviews/` inside; one root `.cmd` each |
