# Gotchas

Traps found in passes 1–110, kept so they are not found again. Read the parts that touch what you are changing. Each
item names the journal entry with the full story where it has one. Add to this file when a new trap is found; remove an
item when the code no longer allows the trap.

## Runs, identity and borrowing

(From the working rules before 26 Nov, unchanged except where marked.)

- **The control is borrowed, not re-run** (user, 21 Nov; journal 135). Runs are deterministic and each sweep's control
  is the previous sweep's adopted arm, so every sweep `.cmd` passes `--control <arm> --identity <the previous adopted
  arm> --borrow-control` to `run_sweep`: the control runs its first seed only (`--check-seeds` for others), compared
  field by field with the reference (`<control>/borrowed.txt`); identical, its other seeds are hard-linked from the
  reference (no extra disk); different, the sweep runs them itself. The comparison still uses six control seeds, and
  the vetoes and every review read the control as an ordinary arm. Hard links share space: deleting an old reference
  sweep frees nothing while a later sweep links its files. The borrow needs the reference run with numerically
  identical code at the sweeps' `threads` 2 (pass 102: a speed-up checked at `threads` 1 broke the borrow at the 14th
  digit; journal 137 §1). Until pass 104 a run at `threads` 2 did not repeat itself (ducc0's threaded adjoint, the slab
  parcels' coefficients from step 7; journal 140 §1), so no borrow had ever succeeded; the adjoint now runs on one
  thread. Check a change claimed "identical" at `threads` 1 and at `threads` 2, and a run against itself at `threads` 2,
  with `precheck-p104-repeat/rep.py` in the device shell (n = 48, 70 steps, two runs side by side, ~2.5 min); the
  sweep's borrow check seed tests it again on Windows.
- **Cloud pre-checks run on seeds outside the sweep's six** (user, 31 Oct): 1007 first, then 1008, ... (1033 is the next unused; 1032 was the pass-110 full-length pre-check, 1031 the pass-109 one, 1030 pass 108's, 1029 pass 107's, 1028 the pass-105 trace and its warm state, re-used in passes 106–108; 2001–2006 are reserved for confirmation) (never used
  before 31 Oct). A cloud run of seed 1001 shares its start with Windows' 1001, so it is a correlated sample; an
  unused seed is an independent one, and the pre-checks stop tuning toward one seed. Run its control beside it (the
  old cloud-1001 numbers do not carry over). Report it beside the six-seed means as a seventh read; it is not pooled
  into a rule fixed on six seeds. Paired probes may use any seed (they compare arms on the same states).
- Runs drift apart from the first steps. Measure an effect smaller than that drift paired, on the same states.
- Step-by-step identity checks need `params={'threads': 1}` (ducc0 threads differ at 1e-12). When editing the default
  path, check that it stays identical step by step (n = 48, ~70 steps). Call `sim.state()` once before the first step
  and after every step, as `simulate_plates` does: `state()` sets `z_last`, which the measured dynamic topography
  reads at the next step (journal 88), so harnesses that call it at another cadence follow another trajectory.
- Sweep folders are named only by the parameters that differ (`MAX_NAME` 80, else `combo-NN`). Shared parameters go in
  `sweep.json`, per-combination ones in `<combo>/params.json`. Runs save `params_effective`. Reviews read older runs
  with `crust_review.LEGACY_DEFAULTS`.
- `planet_meso_check --landscape latest` prefers the newest theta calibration, so pass explicit landscape folders.
- Continental crust in the reviews is label-free: total (crust + sediment) ≥ 20 km on both sides (`crust_review` at
  46 km; `margin_width` on CRUST1.0's 1° cells via `to_crust1_cells`). The flag reading is printed beside it.
- Sweeps save every 5 Myr over the last 40 (`sweep_plates --save-every`, default 5; pass 48). The landscape stage needs
  every Myr: rerun the chosen seed with `simulate_plates` (default 1).
- `PARAMS` are the working baseline (since 29 Oct): a bare run is the current model. New code reads `data/` and
  `outputs/` through `planet_grid.paths` (`DATA`, `OUTPUTS`), never from paths relative to a file.
- Repeat check (pass 104): `precheck-p104-repeat/rep.py <n> <seed> <threads> <steps> <out.json>` writes the md5 of every
  `state()` array per step; run it twice (or on two code versions) and compare the files.
- Pre-check seeds used so far: 1007–1032 (1032 the pass-110 full-length pre-check, 1031 pass 109's, 1030 pass 108's,
  1029 pass 107's, 1028 the pass-105 trace and its warm state). 2001–2006 are confirmation seeds only. Earlier pre-checks
  ran in the cloud workspace, whose trajectories differ from this machine's: their numbers do not carry over.
- **V1 reads the default path only.** The sweep's pytest runs with today's PARAMS, so it does not test the arm. Before
  adopting an arm, run the planet tests with the arm's settings put in PARAMS (pass 111: the Rift Area Sweep's C
  passed every veto as run, then failed 11 tests as the default, two of them an unbooked thinning of oceanic crust).
- Launching a `.cmd` from Claude Code's shells: the tool shell sets `NoDefaultCurrentDirectoryInExePath=1`, so
  `call "Name.cmd"` fails with "is not recognized" even in the right folder; write `call ".\Name.cmd"` (pass 111).

## Watch-points and known defects (STATE §6 as of pass 110, unchanged)

- **The thermal age is a bad state variable where Lamb's L is flat** (pass 111, journal 150). `lith_frame` 'lamb'
  holds L at 114 km over 270–770 Myr and at 178 km from 1,750 Myr (slopes 1e-6 km/Myr), so
  `thickened_therm` → `age_of_thickness` turns any collision gain into a jump to the stretch's end (300 → 770 Myr,
  1,800 → 3,600 Myr for a 0.01% thickening). L barely moves at the jump, but the node then ages along the steeper
  stretch after it, and any averaging of ages (the 'overlap' read) mixes in age space, where L is concave. Under
  'overlap' the continents' mean lithosphere age rises 5.5 Myr per Myr of run (median 3,600 Myr by 0 Myr); under
  'nearest' 1.3. Read the age's clock bound (the mean may rise at most ~1 Myr per Myr without resets) before reading
  L or the datum level. **The effective age is an inverse-table clock, not a measure of support** (review of 27 Nov,
  journal 151 §2): 4 m of thickening on a 40 km column takes 1,800 Myr to 3,600 Myr while L changes by ~0.002 km, and
  material selection moves an area mean. Measure L, column mass and uplift in km and m on persistent parcels.
- **A floor enforced on the averaged input is not enforced on each recipient** (review of 27 Nov; journal 151 §3).
  Under `remap_rule` 'overlap' the return spreads a world cell's change over every plate cell under it. A forearc cut
  limited to the world cell's removable excess then lands on contributors with nothing to remove, and the per-plate
  net correction moves basement that no process touched. Seed 1001 over 250 Myr: 55% / 45% of 29,591 below-floor
  crossings; the correction's share grows late (66% in the last 50 Myr). A global clamp hides this and creates crust
  while passing the floor assertion. Debit only available material, book any shortfall, and test on a small known
  example (`outputs/review-astra/2026-11-27/remap_trace.py` counts the crossings by term).

- **The orogeny-age hand-off** (journal 141 §1 item 4; **fixed in pass 105**, journal 142 §1): NaN in
  `orogeny_age_myr` is the validity mask (no orogeny recorded, or no plate; old crust), a finite 0 an active belt;
  readers must not turn NaN into 0. Every landscape and meso result before 23 Nov carries the old defect (active belts
  as old rock); `--age-rule legacy` reproduces it for replays. `run_planet_landscape` computes the erodibility every
  run (never from a cached forcing) and refuses a forcing folder built for another plate run.
- **The return transfer and the composition's read** (journals 143–144): `composition._update_rasters` ('read') gives a
  step's increments to the plate nodes the world read (the mean when two read one); ~6.8% of continental plate nodes
  get nothing per step, but that is not what roughens the next composition (pass 107: giving them a neighbour's
  increment, the world's end state, or halving the shared reads changes nothing; the forward map is worse because it
  moves ~7% of the world's changes onto a neighbour's material). The nearest-node maps between the world grid and a
  rotated copy cannot be one to one: the cube's cells differ in size across a face, so the rotated grid is denser in
  some places and sparser in others. Every step's surface is a repaired one (erosion takes back ~2/3 of what the
  composition adds).
- **`remap_rule`** (pass 108, off; 'overlap' failed V3 in the Remap Sweep, now in the Rift Area Sweep with `rift_area` 'plates'): the composition's view is an area mean of the plate
  cells under each world cell (crust type and orogeny age stay the nearest node's), so a plate raster read by nearest
  node is rougher than the view (663 → 709 m crust roughness over 20 traced steps under 'overlap', 673 → 772 under
  'tent'); read the plate rasters through the composition, not by nearest node. The return is exact per plate and crust
  type (net, by an additive correction; pass 109: separate gain / loss rescaling blew up to inf), plate crust clamped
  at zero (`return_thick_clamped_m3`), so the crust budget's update part no longer carries nearest's writeback loss
  (−0.146 of crystalline erosion on seed 1031, control +0.011; −0.123 … −0.177 on the Remap Sweep): the rift's opened
  area showed through (pass 110; `rift_area` below). `trace_shore.py --floor`
  (`precheck-p108-floor/`) reads `ret_*` by nearest node, so under a remap only `ret_next` (built by `_world_from`) is
  the composition's read.
- **`rift_area`** (pass 110, off; 'plates' in the Rift Area Sweep): under 'world' the pure-shear zone gives the gap's
  world cell volume while the plates switch on, over the gap, the owning plate's unmasked nodes whose forward nearest
  node is it (~0.91 of the area; 13% of gaps open nothing, about half more than their cell): the update part carries
  the difference (−0.47 km³/yr on the warm state). 'plates' strains by that opened area (`_opened_area`, computed in
  `_rift` with this step's rotations and ownership, the set `_update_rasters` switches on); a gap that opens nothing
  takes no volume and is shown as the zone's crust (Tbar), carried by no plate (`rift_unopened_nodes`). The rift's area
  ledger and `rift_stretched_m3` count the opened area under 'plates'. `rift_pure_shear` takes `opened` (default the
  cell areas: the old kernel bit for bit).
- **A non-finite run stops a sweep at scoring** (`sweep_plates.score`, pass 109): read `run.log` for `nan` and every
  pre-check's log for `RuntimeWarning` (pass 108's two overflow warnings per overlap run were the early sign).
- **The trace** (`precheck-p105-trace/trace_shore.py`): it wraps `_world_from`, `_erode` and `_update_rasters` on the
  instance and restores the sea level after every extra solve (`quiet_state` for `state()` too): traced = untraced
  field by field (`check-identity`). A harness that calls `state()` an extra time moves the trajectory (the next
  solve's loads read the sea level).
- **A one-arm sweep is named `combo-01`** (journal 141 §5): with a single combination nothing differs, so
  `sweep_plates` falls back to `combo-NN`; a `.cmd` that names the control by its grid value (`rift_traction-True`)
  then matches nothing, and the coast reviews, `identity` and any step that names the folder fail. Name the control
  `combo-01` there, or fix the naming before the next one-arm sweep.
- **Selection on elevation** (journal 141 §1): a band picked by low elevation (the drowned normal band) sits on the
  field's lows under any placement (rotations give 4.02–5.70%); read entrants and exits with water re-solved before
  reading a field's share of a band.
- **`rift_traction`** (default since pass 104): the flow's traction enters `_cut_forces` only (the ranking, the
  tension test, the cut chosen), not the split rate (Earth's production / L) or the plate motions (`force_law`
  'lengths'); the velocity is D's scaled per degree like its topography (`MantleFlow.placed_velocity`; the plate-blind
  part of the field has no flow); moot without `dyn_rule` 'slabs'. `PREDATING` False for older runs.
- **Repeatability at `threads` 2** (pass 104): ducc0's threaded adjoint did not repeat run to run (the slab record's
  coefficients from step 7), so no borrowed control had ever borrowed; the adjoint now runs on one thread (`point_alm`,
  `SphereFlexure.deflection`). Check a change on the default path, and a run against itself, at `threads` 1 and 2 in
  the device shell (`precheck-p104-repeat/rep.py`); `Repeatable Control.cmd` tests it on Windows. A new call of a
  ducc0 adjoint (`adjoint_synthesis_general`) takes `nthreads=1`.
- **`flow_placed` 'earth_corr'** (default since pass 102; moot unless `dyn_rule` 'slabs' with `flow_power` 'measured'):
  the plate-blind rest (0.81 of the power, nearly all of l ≥ 3) is selected by the continents like any field (journal
  137 §3). Tests that exercise the hand-over as built pin 'all'. The snapshot's `flow_topography_alm_m` is
  the placed part only; the reviews read `dynamic_topography_air_m` (the whole field with swells) for the spectrum and
  P-b. `EARTH_SLAB_CORR` is noisy at l ≥ 4 (inside the band of unrelated fields; r² ≤ 0.16 there); negatives clip to 0.
  The rest's knots come from the measured field's stream after the start's knots, so any change to the start's draws
  changes the rest too.
- **`flow_power` 'measured'** (default since pass 100, with `dyn_rule` 'slabs'): the record starts empty, so the start's
  field is 'measured' (same stream) and hands over in one `dt_knot_myr`; the review's D bands P-a / P-e do not apply
  (its share is Earth's by construction, its velocity the unscaled flow's); P-c reads nan at the sweeps' save cadence.
  The rescale amplifies the flow's weak high degrees (its P30/P2 ~10⁻³ against Earth's 0.035): the slab pattern stands
  at l ≳ 15, where Earth's signal is mostly lithosphere and upper mantle (Davies et al. 2019). The field puts the
  ridges on highs (+229 m on floor 0–10 Myr; the floor law is Earth's normal floor, which already holds Earth's mean
  residual there) and the continents' interiors on highs (pass 100). Tests that exercise D as built name
  `flow_power` 'flow'.
- **Borrowed control** (since journal 135; never succeeded before pass 104's fix): a borrowed seed's snapshots and `stats.json` are hard links to the
  reference sweep's files; a tool must never rewrite a run's existing `.npz` / `.json` in place (both sweeps would
  change). `<control>/borrowed.json` says which seeds were borrowed; the `identity` review's verdict on them is by
  construction, the check seed's (`borrowed.txt`) is the test. The reference must have run with code that is
  identical at `threads` 2: pass 102's borrow failed at the 14th digit because the reference predated a67e919 (checked
  at `threads` 1 only; journal 137 §1). A code change claimed "identical" is checked at `threads` 2 too.
- **`start_cover` 'accommodation' needs a field at the start**: with a bare level nearly collinear with the crust (no
  dynamic field) its regression weights grow large and Earth's thickest cover lands on thick crust (pass 99's first
  build: a 9.5 km start node). Any rule that changes the start's field must keep one.
- **`marine_diffusion` 'none'** (default since pass 98): the plate layer only (the landscape stage keeps `Kd_marine`);
  deep submerged nodes away from the shelf law and the gravity paths change only by deposits; the gravity paths
  deposit only, so steep submarine slopes are not relaxed (the margins band failed on the means, 201 km / 13%).
- **`start_cover`** (default since pass 97): the bare level is the model's local terms (Airy line before the rim, C, the
  air-loaded field at the start), not the start's isostatic solve with flexure (ranks only); the regression weights
  use the start's own r(crust, bare) (Earth's is +0.405, the start's ~+0.6–0.7); R² capped at 0.99; moot on a bare start.
- **`sediment_rock`** (default since pass 96): the loose cover goes to the kernel in solid metres (the 'cover' law's
  approximation); rock lost per substep by cuts below the loose cover; merges and piece transfers stack rock with the
  sediment; the world's rock is clipped to its sediment. Needs `sediment_entrainment` 'cover' (tests that pin
  'bedrock' pin `sediment_rock` 'none').
- **The sea follows the continents** (pass 95): read land with the sea on the datum and `basin` (young floor).
- **V3 per seed**: a seed's crust residual can leave the control's range from merge steps, holes or the raster update
  alone: read `v3_split` before reading it as a leak (passes 97 and 98 passed by the update-share clause).
- **`backarc_extension`** (default since pass 93): ~20% of the expected area is taken; read `backarc_area_m2` /
  `backarc_expected_m2` beside the thin share; continental upper plates only; its thinning and strips cancel in
  `crust_budget` (the residual shows their jitter); Earth's regime shares rest on regional signs (journal 119 §2).
- **Tests pin what they were set on**: the old margins (`test_lithosphere`'s gate test, `test_torque_law_keeps_the_ledgers`:
  `backarc_extension` False, rollback 17), Artemieva's frame, a bare start (`start_sediment` 'none', journal 116 §7),
  the copula start (`test_start_sediment`'s cover test); `test_isostatic_elevation` reads the run's floor law; the
  Airy-column and weld tests pin `lithosphere` 'none'; exact compares need `threads` 1; tests counting the landscape
  kernel patch the kernel their entrainment rule calls.
- A rule that accumulates state gets a full-length cloud pre-check (journal 102); start it detached (`setsid nohup` from
  a launcher script) when two runs need more than the cloud's 2-hour background limit. The cloud's pytest needs
  `outputs/pytest-tmp` to exist (the conftest's basetemp); without the Earth reference files 11 tests skip.
- Towers: `max_elevation.txt` (V5) on every sweep; split one with `precheck-p65-towers/tower_support.py` (ETOPO's
  highest node 5,747 m). Edge measures at 0 Myr are one snapshot: read the window (`precheck-p71-window/edge_window.py`).
- A new term in `surface.elevation` goes into `elevation_check` in the same pass (journal 92). The young-ocean gate is
  local Airy while the model has flexure (Te 30 km): loaded margins can push a 0 Myr reading past ±300 m (journal 110 §5).
- Reviews read a run's own parameters; older runs take `crust_review.PREDATING` (incl. `start_cover` 'copula',
  `marine_diffusion` 'uniform'): add a key when a PARAMS default changes.
- `isostasy.EARTH_FLOOR_DEPTH_M` depends on the ocean's compaction law and densities: recompute if they change.
- Reading Earth: seismic L in subduction zones may include the slab; CRUST1.0 bins carry its errors (read beside distance
  inland and Lamb's means); the label-free continent includes plateaus and arcs (read small pieces by type); GlobSed's
  medians are flat at ~100 m on 25–95 Myr floor (read means); the thick tail counts crust + sediment; Holdt's field on
  the continents combines 10,874 oceanic and 3,777 continental points in one smooth global fit (journal 141; only seven
  seismic points fall on the 46 km drowned-normal mask), so Earth's z − dyn is a gridded estimate, not an observed
  field-free Earth; `shelf_profile`'s Earth sample leaves ice-sheet nodes out but counts land under ice as land.
- Reading the planet: a plate snapshot's `elevation_m` is relative to sea level; nearest-node sampling loses ~17% of
  continental sediment increments (~4% of pelagic); a node can keep a dropped plate's id (reviews skip it); rifted
  margins read with the flag (`rift_production`); the survival record's 'lost' end includes landings off the crust;
  erosion is measured by recording it (`precheck-p59-export/erosion_record_run.py`); the plates veto is T1 alone; the
  snapshot field `tectonic_thickening_m` is the last step's, not cumulative; `dyn_tail` reads crystalline ≥ 30 km,
  `continent.py` total ≥ 30 km with the flag (their quantiles differ).
- `crust_flow_rule` other than 'fitted': steady geotherm; read the ceiling's share (`crust_flow_capped_edges`).
- **Inland distance depends on what counts as water** (review of 27 Nov). `geography`'s `far300` measures from any
  water node, so isolated depressions make coasts. `far300_sea` measures from the largest connected water body only.
  On A, the low share among far land is 3.6 against Earth's 19.4 on the first and 11.9 against 22.0 on the second.
  Read both side by side. Neither says which depressions should be lakes, dry basins or seas, and the sea solve
  (`isostasy`) counts every node below the level regardless of access.
- **C's coast has costs the crossing count does not show**: shallow water 0–200 m and the shore step (`shelf200`,
  `coast_step` in `geography`). Count them in any terrain claim for C.
- **Elevation at equal crust thickness** (review of 27 Nov): A has more 60–65 km crust than Earth on 1° cells but
  stands 1.2 km lower there, while its 45–50 km crust is already higher than Earth's. A datum shift or general
  thickening moves the wrong part of the distribution. Compare at matched thickness on common 1° support, since
  CRUST1.0 gives one column per 1° cell.
- **Medians at 0 Myr and means over a window answer different questions.** The reviews report six-seed medians at
  0 Myr; STATE and `geography` report six-seed means over −40 … 0 Myr. Don't compare the two directly.

## Later stages (recorded 27 Nov, journal 151 §6; not worked on while the plates are)

- The landscape removes most of the high ground the forcing supports (seed 1001: places above 3 km on the parent go to
  4,448 m with forcing alone and 638 m finished). Judge it on matched locations, not on a percentile: the percentile
  mixes in where the high ground went.
- The landscape applies the *present* erodibility map through its whole 30 Myr window (`scripts/run_planet_landscape.py`).
  The zero-age contract is right; applying today's activity retroactively is a separate assumption.
- On Earth the archived one-step law overpredicts cosmogenic basin denudation 11.4× (median; 8.4× without the slope
  cap). Different averaging times and repeated regional samples: this tests the initial tendency, not a steady-state
  law.
- Meso: Earth-input and planet-input power ratios do not split the planet's gap between causes (journal 146's "about
  half" corrected).

## Queued corrections (STATE §5 as of pass 110, unchanged)

From the independent review; do them with the next audit update or when the tool is next touched:
- The audit's labels. The 3.3× becomes raw / normalized / applied at the same degrees and loading. 'earth_corr''s share
  is a clipped squared correlation with this solver's reconstruction (model-dependent). The normal-ocean age curve is
  regressed on ETOPO minus the model's own terms. Rollback 9 mm/yr is a frame-specific mean, chosen by the regime
  split. GSRM's 290 km half-width is used as a pinning distance. `rift_traction` is normalized by topographic power.
  Holdt 2022 uses 3,777 continental points.
- `cont_spectrum.py` and `earth_cont_corr.py` add `sedthk` to `crsthk` (CRUST1.0's `crsthk` is already total).
- `planet_cont_corr.py`'s "se over runs" treats one trajectory's epochs as independent.
- `three_factor.py` apportions without a water re-solve (descriptive only).
- ETOPO: confirm the local `surface` raster's product identity (it reads sub-ice bed). The Earth JSON's land
  percentiles are unweighted.
- `plate_texture` qualifies windows at ≥ 98% land: report the qualifying area with every texture score.
- The coast reviews and the repeat step of a one-arm sweep (§6).
- From the 27 Nov review, for the audit: `crust_density_rule` 'pressure' treats densification as equilibrium through
  the whole lower-crust share (an upper bound); label it so beside `eclogite_rule` 'foundering' (65 km, τ 10 Myr).

Small items (older): why the thick tail and export fall with the marine term off (pass 98, not read); `crust_budget`
could book the series' crust ledger as its own terms; the start's craton age draws by seed; pelagic rain uniform; D's
slab mass reads GDH1; freeze-run stripes and endorheic basins; the torque law's forces 1.27× too large; a split plate
can run 3 cm/yr above the ceiling for ~10 Myr; the shelf law's κ is band-picked; the graded shelf is assumed; total
sediment 507 Mkm³ against Earth's 512 (the far floor short, the drowned continent over). Run speed (journal 133 §4):
the per-plate nearest-node lookups and running the tests alongside the runs are left as options.
