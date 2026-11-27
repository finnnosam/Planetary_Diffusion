@echo off
rem Rift Speed Sweep (7 Nov; journal 93, seventieth pass). Rule fixed 7 Nov before the runs.
rem  The choosing reading (precheck-p70-edge, on the C arm of sweep-20261107-144712, 0 Myr, K1's instrument): of the
rem  non-passive edge's excess over Earth (165k against 131k km), the transform-near part is mostly the instrument:
rem  Earth's PB2002 steps retyped by the planet's own rule give 42.1k km transform-near (Bird's classes 32.7k), the
rem  planet 47.3k (1.12x). The ridge-near part (41.6k against 24.4k) is unbroken rifts: 31.8k km of edge beside a ridge
rem  with no ocean floor of its own plate within 300 km, against Earth's 14.8k beside continental rifts (CRB); the
rem  ocean-ridge-near edge is Earth's (9.8k against OSR 9.6k). Those rifts were split 8-16 Myr ago (most), open at
rem  11-25 mm/yr (full rate; Earth's near-edge ridges 14) and keep ~10-13 km crust at the axis: pure shear breaks a
rem  rift up after ~3 times the zone half-width (100 km) of opening, which at ~19 mm/yr takes ~11 Myr with an edge.
rem  Earth (Brune et al. 2016, Nature 536, 201): rifts open slowly (< 10 mm/yr) and then accelerate abruptly to a fast
rem  phase (South Atlantic > 35 mm/yr within 6 Myr), breaking up ~10 Myr later. K times its splits by Earth's floor
rem  production beside continent, i.e. the fast phase; the planet's split opens at split_speed_cm_yr 2 (assumed).
rem  This sweep: rift_speed_cm_yr, a split across a continent opens at Earth's fast-phase rate (3.5 cm/yr full);
rem  ocean-only splits keep 2; the halves then relax to the force law as before. No fitted number.
rem  Expected (written before the runs, not a band): unbroken-rift edge 31.8k -> ~17k km, K1 up by ~10-15k km.
rem  Arms (2 x 6 seeds = 12 runs, about 35-45 min); folders rift_speed_cm_yr-...:
rem   0     today's PARAMS (control; must equal lithosphere-thermal_age of sweep-20261107-144712: identity.txt)
rem   3.5   + the fast rift phase
rem  Pre-check (cloud seed 1007, 250 Myr, read only, not pooled; precheck-p70-edge): unbroken-rift edge 43.5k -> 38.2k km
rem   (d < 100 km 19.1k -> 9.2k; axis column p50 12.6 -> 9.2 km); rift pairs split 25-50 Myr ago 3.6k -> 11.3k (the
rem   split's speed relaxes to the force law in ~10 Myr); K1 107.5k -> 102.3k (the trajectories part: plates 33 / 37,
rem   subduction-near +5.7k); highest node 6,858 / 6,663 m; volume -1.3 / -0.4%; land p50 1,129 / 1,236 m. Elevation
rem   gates not readable in the cloud. Smaller than expected on one seed; the six decide.
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds (elevation_check
rem   reads C's term since pass 69); V3 ledgers close (tests); V4 plate_pieces_review T1 detached area <= 0.01 on every
rem   seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: rift_speed_cm_yr 3.5 if its arm passes every veto. Its case is Earth's measured fast rift phase
rem   replacing an assumed 2 cm/yr, not a score.
rem  Scored, reported (no veto): K1 and the edge split (edge_split), the edge's nearest-boundary cause (edge_cause:
rem   ridge / transform / stalled / oblique, Earth by the same rule), the unbroken-rift edge (ridge_near2; Earth CRB
rem   14.8k) and its rift pairs (rift_pairs), K2-K4 (read after by Claude), margins band, rifted production, the rift
rem   and survival records, land p25 / p50 / p90, coast_budget, export, volume, label-free area, T2-T5, plates.
rem  Runs with scripts.run_sweep: tests, the 12 runs and the reviews under one progress line; then the four edge
rem  probes below, each into <sweep>\reviews\*.txt.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P70=%CAL%\precheck-p70-edge
%PY% -u -m scripts.run_sweep --name rift-speed --grid rift_speed_cm_yr=0,3.5 --seeds 1001,1002,1003,1004,1005,1006 --control rift_speed_cm_yr-0 --identity %CAL%\sweep-20261107-144712\lithosphere-thermal_age || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Edge probes on %SW%
%PY% %CAL%\precheck-p69-gate\edge_split.py "%SW%\rift_speed_cm_yr-0" "%SW%\rift_speed_cm_yr-3.5" > "%SW%\reviews\edge_split.txt" 2>&1 || goto fail
%PY% %P70%\edge_cause.py "%SW%\rift_speed_cm_yr-0" "%SW%\rift_speed_cm_yr-3.5" > "%SW%\reviews\edge_cause.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%SW%\rift_speed_cm_yr-0" > "%SW%\reviews\ridge_near2_0.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%SW%\rift_speed_cm_yr-3.5" > "%SW%\reviews\ridge_near2_3.5.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%SW%\rift_speed_cm_yr-0" > "%SW%\reviews\rift_pairs_0.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%SW%\rift_speed_cm_yr-3.5" > "%SW%\reviews\rift_pairs_3.5.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
