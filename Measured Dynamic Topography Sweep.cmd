@echo off
rem Measured Dynamic Topography Sweep (1 Nov; journal 80, fifty-eighth pass). Audit step 5 / STATE next item 1: the
rem  dynamic topography from Earth's measured spectrum replaces the fitted random field (dyn_std_m 500, flat power to
rem  degree 30). dyn_rule 'measured' (planet_tectonics.mantle): plate-blind knot fields every 40 Myr with Valentine &
rem  Davies (2020)'s residual-topography power per degree l 1-30 as they publish it (expected posterior power of their
rem  coefficients for the high-accuracy oceanic spot data: rms 0.58 km water-loaded, 0.41 km air-loaded, 0.29 of the
rem  power at l 1-3; their abstract's "~85% in degrees 1-3" is the share of degrees 1-15, journal 80 section 1), held
rem  air-loaded and water-loaded under the sea as D is, blended at constant variance (no pulse). No fitted number.
rem  Arms (4 x 6 seeds = 24 runs, about 50-70 min); folders combo-01..04 (params.json):
rem   combo-01 crust_kappa_m2_yr 0, trench_polarity pair, trench_flux consumed, arc_H_ocean_m 1000, arc_accretion True,
rem       lithosphere none, dyn_rule random: combo-01 of sweep-20261101-165624 rerun; must equal it (identity.txt)
rem   combo-02 combo-01 with dyn_rule measured
rem   combo-03 combo-02 with lithosphere thermal_age (C read again)
rem   combo-04 today's PARAMS with dyn_rule measured (control: combo-01 of sweep-20261101-115244, not rerun)
rem  Rule (fixed 1 Nov before the runs), six-seed means unless stated:
rem  Part 1, the measured field (combo-04 for the default; combo-02 reported the same way):
rem   build checks (dynamic_topography, Earth's observed spectrum air-loaded): M-a share at l 1-3 0.14-0.43; M-b rms
rem   324-507 m; M-c P30/P2 <= 0.1; M-d no 10-Myr rms step > 20%; M-e rate p50 <= 20 m/Myr.
rem   integrity: pytest passes; elevation gates +-300 m on 6 of 6 seeds (elevation_check); plates at 0 Myr <= 40 on
rem   every seed; ledgers close (tests).
rem   If the build checks and the integrity hold on combo-04: PARAMS dyn_rule 'measured' is adopted. It replaces a
rem   fitted value with a measurement, so no terrain band is asked of it (realism rule: no fitted value is kept for a
rem   band). Export 4-8 km3/yr, volume >= -14%, margins and the hypsometry are read against the control and reported;
rem   a guard that holds in the control and fails here is what the fitted field was carrying, and is named as such.
rem   Build checks fail -> the build is wrong (fix, not adopted). A gate fails -> not adopted; where the field moves
rem   that gate is the next step.
rem  Part 2, C read again (combo-03 against combo-02): adopted (PARAMS lithosphere thermal_age, crust_kappa_m2_yr 0 and
rem   B-d2 arm C's settings) only if Part 1 adopts the field and: C1 (coast_budget, Earth on full-land 1-degree cells:
rem   -524 m, band -262 ... -786 m) passes; C2 rise of land inland (coast_profile, <= 772 m) passes; >= 3 of P1-P4
rem   (as the Lithosphere Age Sweep: P1 > 45 km <= 0.13, P2 > 55 km <= 0.036, P3 slope of land > 2 km p50 7.0-17.6
rem   m/km, P4 > 45 km depth p50 58-146 km); C's term guard (coast to interior -119 ... -357 m); and every guard:
rem   gates 6 of 6, volume >= -14%, export 4-8 km3/yr, margins pass, plates <= 40.
rem   Not adopted -> the next step in every branch is export (the narrow belts' erosion rate, STATE item 2); C1 sets
rem   its base: C1 passes on combo-03 -> combo-03; C1 fails -> combo-02, with the interior's remaining excess carried as
rem   the residual (journal 79 section 5).
rem  Pre-check (cloud seed 1008, combo-01 against combo-02; journal 80 section 4) is reported as a seventh read.
rem  Reported, not scored: land p25 / p50 / p90, land fraction, continental crust drowned, area, sea level on the datum,
rem  lowland beyond 300 km, C1 without the dynamic topography and the field's interior lift (coast_budget), T2-T5
rem  (plate_pieces_review), arc growth A1, the applied field's rms.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py tests\test_lithosphere.py tests\test_dyn_measured.py > %CAL%\measured-dyn-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\measured-dyn-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid crust_kappa_m2_yr,trench_polarity,trench_flux,arc_H_ocean_m,arc_accretion,lithosphere,dyn_rule=0:pair:consumed:1000:True:none:random,0:pair:consumed:1000:True:none:measured,0:pair:consumed:1000:True:thermal_age:measured,1000:node:band:300:False:none:measured || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\measured-dyn-pytest.txt "%SW%\reviews\pytest.txt" >nul
%PY% tools\compare_runs.py "%SW%\combo-01" %CAL%\sweep-20261101-165624\combo-01 > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% -u -m planet_reviews.dynamic_topography "%SW%" %CAL%\sweep-20261101-115244\combo-01 --json "%SW%\reviews\dynamic_topography.json" > "%SW%\reviews\dynamic_topography.txt" 2>&1
type "%SW%\reviews\dynamic_topography.txt"
%PY% -u -m planet_reviews.coast_budget "%SW%" %CAL%\sweep-20261101-115244\combo-01 --json "%SW%\reviews\coast_budget.json" > "%SW%\reviews\coast_budget.txt" 2>&1
type "%SW%\reviews\coast_budget.txt"
%PY% -u -m planet_reviews.coast_profile "%SW%" %CAL%\sweep-20261101-115244\combo-01 --json "%SW%\reviews\coast_profile.json" > "%SW%\reviews\coast_profile.txt" 2>&1
type "%SW%\reviews\coast_profile.txt"
%PY% -u -m scripts.plate_pieces_review --latest --control combo-01 > "%SW%\reviews\plate_pieces_review.txt" 2>&1
type "%SW%\reviews\plate_pieces_review.txt"
for %%r in (crust_review crust_budget elevation_check margin_width interior_probe) do (
  %PY% -u -m scripts.%%r --latest > "%SW%\reviews\%%r.txt" 2>&1
  type "%SW%\reviews\%%r.txt"
)
echo Reviews saved in %SW%\reviews
echo done
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
