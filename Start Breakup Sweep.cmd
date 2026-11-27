@echo off
rem Start Breakup Sweep (6 Nov; journal 88, sixty-sixth pass). Rule fixed 6 Nov before the runs.
rem  Read of the Breakup History Sweep (sweep-20261105-211423): V1-V4 pass on every arm; V5 (towers) passes on combo-02
rem  and combo-05 only (combo-03 10,857 m on seed 1003, combo-04 9,705 m on 1001, combo-06 on 1001 / 1002 / 1006). By
rem  the rule's order combo-05 is adopted: collision_spread earth, eclogite_rule foundering, shortening_match True,
rem  split_rule earth, dyn_rule measured are PARAMS now. Scored: K makes 1.75x Earth's rifted margin per unit continent
rem  (control 1.14x) and still none of Earth's quiet 250-200 Ma: every arm makes 4,150 / 4,210 km/Myr at -240 / -230 Myr
rem  (Earth 977 / 399), because the start lies on five continental plates moving apart (start_breakup 'dispersal'), and
rem  10-12% of the supercontinent on oceanic plates whose Voronoi cells reach past its coast. K1 / K2 / K4 fail on every
rem  arm (passive edge 86-106k km at 1 deg), K3 passes: the passive edge is set by survival.
rem  Built this pass (default off; the new PARAMS path identical step by step to combo-05's settings on the old code,
rem  n = 48, 2 seeds x 70 steps, threads 1, recorder running):
rem   start_breakup whole: the whole supercontinent (and the continental plates' cells) is one plate with a random pole at
rem     2-4 cm/yr, drawn from the same numbers; it breaks only where split_rule earth splits it. No new number.
rem   planet_reviews.passive_survival: every run now records, each Myr, what ends each passive-edge node (carried with
rem     its plate), and the births (output only). Earth: Bradley 2008, mean lifespan 181 Myr (Permian-Neogene 130).
rem  Arms (4 x 6 seeds = 24 runs, about 60-75 min); folders combo-01..04 (params.json):
rem   combo-01 today's PARAMS (control; must equal combo-05 of sweep-20261105-211423: identity.txt)
rem   combo-02 + start_breakup whole
rem   combo-03 + start_breakup whole + start_plates coast (the supercontinent's rim a plate boundary all round; rejected
rem            26 Oct alone on T5 and the gates, journal 56; T5 is scored now)
rem   combo-04 combo-03 + start_poles forces (oceanic start plates move by the force law with the rim as trench: the
rem            circum-Pangea subduction girdle, Collins 2003; Cawood 2005)
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2 elevation gates +-300 m on 6 of 6 seeds (elevation_check);
rem   V3 ledgers close (tests); V4 no runaway: plate_pieces_review T1 detached area <= 0.01 on every seed; V5 no towers:
rem   on every seed the highest node over the run (max_elevation.txt) <= 8,850 m or no higher than combo-01's there.
rem  The arms' own test S1 (what they are for): rifted production (rift_production.txt, six-seed mean) over the bins
rem   -240 ... -200 Myr <= 1,000 km/Myr on average (Earth 977 / 399 / 242 / 0 / 0, mean 324).
rem  Adoption: of the arms passing every veto and S1, combo-04 (the girdle and the forces) before combo-03, then
rem   combo-02. Its settings become PARAMS. None -> nothing adopted; the next build targets what failed, named in the
rem   reading.
rem  Scored, reported (no veto): S2 the decade of peak rifted production inside -160 ... -100 Myr (Earth's 160-110 Ma);
rem   rifted production per unit continent against Earth's (1.0x); survival (passive_survival.txt, last 150 Myr):
rem   lifespan = passive edge / edge ending per Myr against 130-181 Myr, and what ends it; read after by Claude (journal
rem   85 probes): K1 passive edge (1 deg) >= 138k
rem   km, K2 50-200 Myr margins >= 123k km, K3 margin age p50 85-141 Myr, K4 crust 27.5-32.5 km >= 27.0 Mkm2 and drowned
rem   continent >= 100 km2 per km of coast; export, volume, margins band, land p25 / p50 / p90, > 45 / > 55 km,
rem   label-free area, sea on the datum, C1, T2-T5, foundered volume (crust_budget, ~1.1 km3/yr on Earth).
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py tests\test_lithosphere.py tests\test_dyn_measured.py tests\test_margin_wedge.py tests\test_collision_spread.py tests\test_breakup_match.py tests\test_start_breakup.py > %CAL%\start-breakup-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\start-breakup-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid start_breakup,start_plates,start_poles=dispersal:voronoi:random,whole:voronoi:random,whole:coast:random,whole:coast:forces || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\start-breakup-pytest.txt "%SW%\reviews\pytest.txt" >nul
%PY% tools\compare_runs.py "%SW%\combo-01" %CAL%\sweep-20261105-211423\combo-05 > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% tools\max_elevation.py "%SW%" > "%SW%\reviews\max_elevation.txt" 2>&1
type "%SW%\reviews\max_elevation.txt"
%PY% -u -m planet_reviews.rift_production "%SW%" --json "%SW%\reviews\rift_production.json" > "%SW%\reviews\rift_production.txt" 2>&1
type "%SW%\reviews\rift_production.txt"
%PY% -u -m planet_reviews.passive_survival "%SW%" --json "%SW%\reviews\passive_survival.json" > "%SW%\reviews\passive_survival.txt" 2>&1
type "%SW%\reviews\passive_survival.txt"
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
