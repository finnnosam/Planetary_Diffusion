@echo off
rem Craton Lithosphere Sweep (7 Nov; journal 90, sixty-eighth pass). Rule fixed 7 Nov before the runs.
rem  Read of the Initiation Site Sweep (sweep-20261107-094624): V1-V5 pass on every arm (pytest 140; gates 6 of 6; T1
rem  <= 0.0021; highest node <= 8,183 m); the control equals combo-04 of sweep-20261106-211658 (identity.txt). D1:
rem  init_site ocean passes every veto -> adopted (PARAMS init_site 'ocean'). D2: B-d3 fails its own test O1 (ocean-
rem  ocean convergence without a slab 1.3k km against >= 8.7k; Bird 17.4k) -> not adopted. Rift record (six seeds, last
rem  150 Myr): young rifted margin enters 1,940 km/Myr and 0.50 of it becomes passive in ~13 Myr (Earth ~0.6-0.85),
rem  0.19 ends at a trench; the survival record's births stay 428 km/Myr (journal 90 section 3).
rem  This sweep: the land's height (STATE item 1, the scorecard's largest failure: land p50 1,295 m against 384; the
rem  interior's z - Airy +209 m against Earth's -367 m at 600-1,200 km inland). C (lithosphere 'thermal_age', built
rem  1 Nov, journal 78): each continental node's tectonothermal age sets its thermal lithosphere, L = 93.6 + 0.04 t km
rem  (Artemieva 2006), whose excess density 30 kg/m3 (Lamb et al. 2020) stands the surface at
rem  -(30 / 3300)(L - 130.67 km) about the Airy line; L_ref is Earth's continental mean on the same age map, so the
rem  term moves no volume on Earth's age distribution. No fitted number. Not adopted twice (journals 79, 81) by C1/C2
rem  bands on models whose crust was too thick (interior 49-52 km); the crust is now Earth's inland (40 km) and the
rem  rule since 4 Nov adopts a realistic mechanism on the vetoes, scoring the terrain.
rem  Arms (2 x 6 seeds = 12 runs, about 35-45 min); folders lithosphere-...:
rem   none          today's PARAMS (control; must equal init_site-ocean_subduction_initiation-instant of
rem                 sweep-20261107-094624 on the shared values: identity.txt)
rem   thermal_age   + C
rem  Pre-check (cloud seed 1010, 250 Myr, read only, not pooled; precheck-p68-craton): C lowers land p25 / p50 / p90
rem   782 / 1,319 / 3,004 -> 555 / 1,110 / 2,434 m; interior (600-1,200 km) z 1,006 -> 854 m; C's term coast to
rem   interior -152 m (its guard -357 ... -119); C1 -131 -> +85 m (without dyn -21 -> -119); highest node 7,758 -> 6,848 m;
rem   volume -2.4 -> -0.5%. The far interior (> 1,200 km) is 283 Myr old with C against 1,988 in the control's age map:
rem   collisions reset it (one seed; the trajectories part from the start). Elevation gates not readable in the cloud.
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2 elevation gates +-300 m on 6 of 6 seeds (with C the
rem   thin-crust offset turns positive: STATE watch-point); V3 ledgers close (tests); V4 plate_pieces_review T1 detached
rem   area <= 0.01 on every seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the
rem   control's there.
rem  Adoption: C (lithosphere 'thermal_age') if the thermal_age arm passes every veto. Its case is Earth's measured
rem   lithosphere (age-thickness and density), not a score.
rem  Scored, reported (no veto): C1 (z - Airy interior minus coast; with and without the dynamic topography) and the
rem   interior's z - Airy (coast_budget), C's applied term by band, land p25 / p50 / p90, land fraction, continental
rem   crust drowned, shelf, sea on the datum, export, volume, label-free area, > 45 / > 55 km, margins band, K1-K4 (read
rem   after by Claude), the rift and survival records, T2-T5, foundered volume.
rem  Runs with scripts.run_sweep (pass 68): tests, the 12 runs and the reviews under one progress line (elapsed against
rem  predicted time); each step's time goes to outputs\planet-calibration\timings.jsonl for the next prediction.
rem  Printouts: the tests in outputs\planet-calibration\craton-lithosphere-pytest.txt, each run in <arm>\seed<s>\run.log,
rem  the reviews in <sweep>\reviews\*.txt, predicted against actual times in <sweep>\reviews\timings.txt.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
%PY% -u -m scripts.run_sweep --name craton-lithosphere --grid lithosphere=none,thermal_age --seeds 1001,1002,1003,1004,1005,1006 --control lithosphere-none --identity outputs\planet-calibration\sweep-20261107-094624\init_site-ocean_subduction_initiation-instant || goto fail
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
