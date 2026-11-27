@echo off
rem Sediment Cover Sweep (15 Nov; journal 110, eighty-fifth pass). Rule fixed 15 Nov before the runs.
rem  Why: STATE section 0, items 2-3, the sea's supply. The planet erodes Earth's gross volume from its land but keeps
rem  ~90% of it on its continents (export 2.0 against Earth's 4-8 km3/yr), and keeps it high: continental sediment on
rem  land 500-3,000 m 116 against Earth's 41 Mkm3, on the shelves and drowned margins a third of Earth's (Thick Tail
rem  Sweep, crust_flow_rule-rey). The 46 km erosion law's incision threshold theta (200, from OCTOPUS 10Be) holds on
rem  every surface: a node below it only receives sediment, so loose deposits on plains and plateaus never leave by
rem  rivers. On Earth the threshold is bedrock's (the stream power a river needs to cut rock); loose sediment is
rem  entrained and carried on by far weaker flows (alluvial rivers are transport-limited). sediment_entrainment 'cover'
rem  applies the threshold to bedrock only: the sediment layer is entrained by the same law at theta 0 down to its
rem  base (graph_landscape.erode_and_deposit_cover). No new number.
rem  Arms (2 x 6 seeds = 12 runs); folders sediment_entrainment-...:
rem   bedrock  today's PARAMS (control; must equal crust_flow_rule-rey of sweep-20261114-210918: identity.txt)
rem   cover    the threshold on bedrock only
rem  Rule. Vetoes on the cover arm: V1 pytest passes (incl. test_sediment_cover); V2' elevation gates +-300 m on 6 of 6
rem   seeds; V3 ledgers close (crust_budget residuals of the control's size); V4 plate_pieces_review T1 detached area <=
rem   0.01 on every seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: cover if it passes every veto (by mechanism, not by score); else bedrock kept and the failing veto is
rem   the next thing to work on.
rem  Full-length cloud pre-check, seed 1016 (unused), bedrock / cover, read only, not pooled (journal 110 section 5).
rem  Expected (not a band): export and sediment past the shelves up, sediment on land 500-3,000 m down, the sea up
rem   (shallower), land p50 down by the sea's share; continental volume at most a little lower (sediment leaves the
rem   continents for the floor, where trenches take some).
rem  Scored, reported (no veto): crust_budget (export, to the floor, deep-sea fill), sediment_where, sediment_by_height,
rem   three_factor, basin, continent, hyps, column_split per arm, thick_tail, lith_crust2; land fraction and
rem   p25 / p50 / p90, crust bands, volume, K1-K4, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
set P83=%CAL%\precheck-p83-sediment
set P84=%CAL%\precheck-p84-thicktail
%PY% -u -m scripts.run_sweep --name sediment-cover --grid sediment_entrainment=bedrock,cover --seeds 1001,1002,1003,1004,1005,1006 --control sediment_entrainment-bedrock --identity %CAL%\sweep-20261114-210918\crust_flow_rule-rey || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P80%\lith_crust2.py "%SW%\sediment_entrainment-bedrock" "%SW%\sediment_entrainment-cover" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\sediment_entrainment-bedrock" "%SW%\sediment_entrainment-cover" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
for %%a in (bedrock cover) do (
  %PY% %P82%\three_factor.py "%SW%\sediment_entrainment-%%a" > "%SW%\reviews\three_factor_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\basin.py "%SW%\sediment_entrainment-%%a" > "%SW%\reviews\basin_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\continent.py "%SW%\sediment_entrainment-%%a" > "%SW%\reviews\continent_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\hyps.py "%SW%\sediment_entrainment-%%a" > "%SW%\reviews\hyps_%%a.txt" 2>&1 || goto fail
  %PY% %P75%\column_split.py "%SW%\sediment_entrainment-%%a" > "%SW%\reviews\column_split_%%a.txt" 2>&1 || goto fail
)
%PY% %P83%\sediment_where.py "%SW%\sediment_entrainment-bedrock" "%SW%\sediment_entrainment-cover" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\sediment_entrainment-bedrock" "%SW%\sediment_entrainment-cover" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
