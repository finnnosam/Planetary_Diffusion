@echo off
rem Breakup History Sweep (4 Nov; journal 87, sixty-fifth pass). Rule fixed 4 Nov before the runs.
rem  Read of the Collision Spread Sweep (sweep-20261104-214903): every arm passed V1-V4, none V5 (towers); nothing
rem  adopted. The no-flow arms' towers are Airy columns (85-112 km crust over 50-500k km2) with R at Earth's foundering
rem  flux (0.5-1.0 km3/yr against ~1.1): without the collapse flow no measured sink holds thick crust (no-flow closed).
rem  Every arm's peaks also carry 1.6-7.8 km of elastic support under node-scale crust spikes; per-step probes (cloud
rem  seeds 1007-1008) found the feeders: the one-node collision dump (fixed by collision_spread earth) and the margin
rem  shortening, which put a trench segment's whole foreland volume on the few nodes of a missing belt (+20-120 km
rem  in one Myr; an along-strike pool tried first did not remove it, journal 87).
rem  Built this pass (default off; default path identical step by step, n = 48, 2 seeds x 70 steps, threads 1):
rem   shortening_match True: the belt is no longer renormalised to the continental nodes a segment has; each belt node
rem     takes the trench flux at its own Gaussian weight and the segment moves the smaller of what the foreland gives
rem     and the belt takes (volume conserved). A missing belt (oceanic forearc) moves little. No new number.
rem   K, split_rule earth: splits across a continent at Earth's rifted-margin production since Pangea (Mueller et al.
rem     2019 age grids, continent beside floor <= 10 Myr, > 300 km from subduction zones; km per Myr per Mkm2 of
rem     continent, by time since the start mapped onto 250 Ma ... 0). The control already makes Earth's total (1.14x per
rem     area, flag instrument) but front-loads it at -240 Myr; Earth's peak is 160-110 Ma.
rem  Arms (6 x 6 seeds = 36 runs, about 90-110 min); folders combo-01..06 (params.json):
rem   combo-01 today's PARAMS (control; must equal combo-01 of sweep-20261104-214903: identity.txt)
rem   combo-02 + collision_spread earth + eclogite_rule foundering + shortening_match True (the towers' feeders)
rem   combo-03 + split_rule earth (K alone)
rem   combo-04 combo-02 + split_rule earth
rem   combo-05 combo-04 + dyn_rule measured (exploratory, adopted by the same rule)
rem   combo-06 combo-04 + collapse_thick_m 60 km (exploratory: the collapse flow starts where Earth's lower crust is
rem            seen to flow, plateau crust of 60-80 km in Tibet and the Altiplano, not at the Tibet-tuned 45 km; the
rem            no-flow arms scored better on area, > 45 km, volume and sea but built towers)
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2 elevation gates +-300 m on 6 of 6 seeds (elevation_check);
rem   V3 ledgers close (tests); V4 no runaway: plate_pieces_review T1 detached area <= 0.01 on every seed; V5 no towers:
rem   on every seed the highest node over the run (max_elevation.txt) <= 8,850 m or no higher than combo-01's there.
rem  Adoption: of the arms passing every veto, combo-05 (also replaces the fitted random field) before combo-06 (the
rem   flow's onset from observation), combo-04, combo-02 (integrity defects removed with measured rules), then
rem   combo-03. Its settings become PARAMS.
rem   None passes -> nothing adopted; the next build targets the veto that failed, named in the reading.
rem  Scored, reported (no veto): export (4-8 km3/yr), volume, margins band, land p25 / p50 / p90, > 45 / > 55 km,
rem   label-free area, sea on the datum, C1, T2-T5, foundered volume (~1.1 km3/yr), rifted production by time against
rem   Earth (rift_production.txt), and read after by Claude: K1 passive edge (1 deg) >= 138k km, K2 50-200 Myr margins
rem   >= 123k km, K3 margin age p50 85-141 Myr, K4 crust 27.5-32.5 km >= 27.0 Mkm2 and drowned continent >= 100 km2
rem   per km of coast (journal 85 probes); the passive margins' survival (what ends them).
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py tests\test_lithosphere.py tests\test_dyn_measured.py tests\test_margin_wedge.py tests\test_collision_spread.py tests\test_breakup_match.py > %CAL%\breakup-history-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\breakup-history-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid collision_spread,eclogite_rule,shortening_match,split_rule,dyn_rule,collapse_thick_m=node:none:False:poisson:random:45000,earth:foundering:True:poisson:random:45000,node:none:False:earth:random:45000,earth:foundering:True:earth:random:45000,earth:foundering:True:earth:measured:45000,earth:foundering:True:earth:random:60000 || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\breakup-history-pytest.txt "%SW%\reviews\pytest.txt" >nul
%PY% tools\compare_runs.py "%SW%\combo-01" %CAL%\sweep-20261104-214903\combo-01 > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% tools\max_elevation.py "%SW%" > "%SW%\reviews\max_elevation.txt" 2>&1
type "%SW%\reviews\max_elevation.txt"
%PY% -u -m planet_reviews.rift_production "%SW%" --json "%SW%\reviews\rift_production.json" > "%SW%\reviews\rift_production.txt" 2>&1
type "%SW%\reviews\rift_production.txt"
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
