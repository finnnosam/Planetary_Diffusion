@echo off
rem Collision Spread Sweep (4 Nov; journal 86, sixty-fourth pass). First sweep under the working rules of 4 Nov (user):
rem  vetoes are integrity only, terrain measures are scored; a pass ends with something built and swept or adopted;
rem  sweeps run even when an arm is not expected to help (the information is the point).
rem  Built this pass (both default off, default path identical step by step):
rem   collision_spread 'earth': a collision's volume is spread over the receiving plate's continental crust with
rem     Earth's measured kernel (crustal thickening rate by distance from Bird's CCB steps, GSRM v2.1; half within
rem     ~290 km, 90% within ~1,200 km) instead of landing on one node (+10-40 km of crust per Myr).
rem   eclogite_rule 'foundering' (R): crust thicker than 65 km (Chapman et al. 2019) loses its excess with a 10 Myr time
rem     constant (Krystopowicz & Currie 2013). No fitted number in either.
rem  Why: the no-flow bundle (collapse flow off + B-d2 arm C) cannot be adopted as it is: its columns reach 90-115 km
rem   and peaks 12-20 km on every seed (combo-01 of sweep-20261101-191306); one-node collision volume is the cause.
rem  Arms (5 x 6 seeds = 30 runs, about 75-90 min); folders combo-01..05 (params.json):
rem   combo-01 today's PARAMS (control; must equal combo-01 of sweep-20261101-115244: identity.txt)
rem   combo-02 today's PARAMS + collision_spread earth + eclogite_rule foundering
rem   combo-03 the no-flow bundle (crust_kappa_m2_yr 0, trench_polarity pair, trench_flux consumed, arc_H_ocean_m 1000,
rem            arc_accretion True) + collision_spread earth
rem   combo-04 combo-03 + eclogite_rule foundering
rem   combo-05 combo-04 + dyn_rule measured (exploratory: run for the information, adopted by the same rule)
rem  Rule (fixed 4 Nov before the runs). Vetoes, every one on its arm:
rem   V1 pytest passes; V2 elevation gates +-300 m on 6 of 6 seeds (elevation_check); V3 ledgers close (tests);
rem   V4 no runaway: plate_pieces_review T1 detached area <= 0.01 on every seed; V5 no towers: on every seed the highest
rem   node over the run (max_elevation.txt) <= 8,850 m (Everest; ETOPO's highest n = 217 node 5,747 m) or no higher than
rem   combo-01's on that seed (today's model reaches 11.2 km on seed 1004).
rem  Adoption: of the arms passing every veto, the most realistic: combo-05 (two fitted values replaced: collapse flow,
rem   random field) before combo-04, combo-03, then combo-02 (both new rules are measured). Its settings become PARAMS.
rem   None passes -> nothing adopted; the next build targets the veto that failed, named in the reading.
rem  Scored, reported (no veto): export (4-8 km3/yr), volume (>= -14%), margins band, land p25 / p50 / p90, > 45 / > 55 km
rem   shares, label-free area, sea on the datum, C1, T2-T5, foundered volume (Earth ~1.1 km3/yr, Clift et al. 2009),
rem   passive edge and margin age (journal 85 probe, read after), K1-K4 readings.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py tests\test_lithosphere.py tests\test_dyn_measured.py tests\test_margin_wedge.py tests\test_collision_spread.py > %CAL%\collision-spread-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\collision-spread-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid crust_kappa_m2_yr,trench_polarity,trench_flux,arc_H_ocean_m,arc_accretion,collision_spread,eclogite_rule,dyn_rule=1000:node:band:300:False:node:none:random,1000:node:band:300:False:earth:foundering:random,0:pair:consumed:1000:True:earth:none:random,0:pair:consumed:1000:True:earth:foundering:random,0:pair:consumed:1000:True:earth:foundering:measured || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\collision-spread-pytest.txt "%SW%\reviews\pytest.txt" >nul
%PY% tools\compare_runs.py "%SW%\combo-01" %CAL%\sweep-20261101-115244\combo-01 > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% tools\max_elevation.py "%SW%" > "%SW%\reviews\max_elevation.txt" 2>&1
type "%SW%\reviews\max_elevation.txt"
%PY% -u -m planet_reviews.coast_budget "%SW%" "%SW%\combo-01" --json "%SW%\reviews\coast_budget.json" > "%SW%\reviews\coast_budget.txt" 2>&1
type "%SW%\reviews\coast_budget.txt"
%PY% -u -m planet_reviews.coast_profile "%SW%" "%SW%\combo-01" --json "%SW%\reviews\coast_profile.json" > "%SW%\reviews\coast_profile.txt" 2>&1
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
