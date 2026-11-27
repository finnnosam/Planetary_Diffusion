@echo off
rem Collapse Flow Sweep (1 Nov; journal 76, fifty-fourth pass). Item 1 (continental crust distribution).
rem  Pass 53 (journal 75): the continents are domes; the > 45 km crust is built by collision and never removed.
rem  Pass 54 found why on cloud seed 1007 (a diagnostic, not a selection): the collapse flow (crust_kappa_m2_yr 1000 with
rem  the e-fold above collapse_thick_m 45 km; Tibet-tuned, a fitted M-target in the audit) spreads collided crust into
rem  broad plateaus that stall just above 45 km, where slopes are low and erosion is slow. With crust_kappa_m2_yr 0
rem  (the measured cold-crust flow, crust_kappa_normal_m2_yr 30, stays) on seed 1007: > 45 km 0.246 -> 0.114, > 55 km
rem  0.098 -> 0.020, slope of land > 2 km 2.6 -> 8.8 m/km, > 45 km depth from its edge 285 -> 109 km; but crust
rem  volume -11% -> -26% (narrow belts erode 2.7x faster and nothing adds crust back: the arc addition is 0.7 km3/yr
rem  against Earth's 3.8, Clift et al. 2009). Arm C adds B-d2 arm C's volume source (pass 41).
rem  Arms (3 x 6 seeds = 18 runs, about 1.5-2 h), everything else PARAMS; folders combo-01..03 (params.json):
rem   combo-01 crust_kappa_m2_yr 1000, trench_polarity node, trench_flux band, arc_H_ocean_m 300, arc_accretion False:
rem       control (today's model); must equal arm A of sweep-20261031-154044 field by field (identity.txt)
rem   combo-02 crust_kappa_m2_yr 0 (collapse flow off), the rest as the control
rem   combo-03 crust_kappa_m2_yr 0 + trench_polarity pair, trench_flux consumed, arc_H_ocean_m 1000, arc_accretion True
rem  Rule (fixed 1 Nov before the runs), six-seed means unless stated:
rem   P1 > 45 km share of the label-free continent <= 0.13 (crust_review's band; Earth 0.086)
rem   P2 > 55 km share <= 0.036 (2 x Earth's 0.018; crust_review)
rem   P3 slope of land > 2 km p50 in 7.0-17.6 m/km (interior_probe E band; Earth 11.1)
rem   P4 > 45 km depth from its edge p50 in 58-146 km (interior_probe E band; Earth 92)
rem   P5 rise of land inland (median at 600-1200 km minus 0-50 km) <= 2 x Earth's (coast_profile; Earth 386 m)
rem   Guards (all must hold): pytest passes; elevation gates +-300 m on 6 of 6 seeds (elevation_check); continental
rem   volume change >= -14% (crust_budget; the E guard of pass 26); export 4-8 km3/yr (crust_budget); margins pass
rem   (margin_width); plates at 0 Myr <= 40 on every seed; ledgers close (tests).
rem   An arm is adopted if it passes >= 3 of P1-P5 and every guard. If both qualify, combo-03 (it carries the volume
rem   source). If neither: nothing adopted; the failing guard sets the next step (volume source first if it is volume).
rem  Reported, not scored: land p25 / p50 / p90, land fraction, drowned crust, area, sea level on the datum, lowland
rem  share beyond 300 km, T2-T5 (plate_pieces_review), sediment subducted, arc growth A1, crust volume by seed.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py > %CAL%\collapse-flow-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\collapse-flow-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid crust_kappa_m2_yr,trench_polarity,trench_flux,arc_H_ocean_m,arc_accretion=1000:node:band:300:False,0:node:band:300:False,0:pair:consumed:1000:True || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\collapse-flow-pytest.txt "%SW%\reviews\pytest.txt" >nul
%PY% tools\compare_runs.py "%SW%\combo-01" %CAL%\sweep-20261031-154044\dyn_rule-random_force_law-lengths_basal_traction-False > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% -u -m planet_reviews.coast_profile "%SW%" --json "%SW%\reviews\coast_profile.json" > "%SW%\reviews\coast_profile.txt" 2>&1
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
