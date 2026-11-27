@echo off
rem Margin Rise Sweep (13 Nov; journal 108, eighty-third pass). Rule fixed 13 Nov before the runs.
rem  Why: STATE section 0, the sea, step (b): the floor's sediment. Measured this pass on today's model (ocean_depth_rule-
rem  earth of sweep-20261113-201653; precheck-p83-sediment/sediment_where.py, sediment_by_height.py): the planet holds
rem  Earth's sediment on land above the sea (153 against 176 Mkm3 bulk) but a third of Earth's on its drowned continent
rem  (45 against 144) and 0.2-0.6 x Earth's on the floor within 1,000 km of the continent (73 against 166 Mkm3; 200-500 km
rem  out 221 m against 1,055 m: Earth's continental rise). S, the margin wedge (built 3 Nov, journal 83; never swept): what
rem  overflows the graded shelf is deposited in front of its source up to Earth's measured floor seaward of passive edges
rem  (ETOPO p50 by distance from the label-free continent, > 300 km from boundaries), before any rest ponds in the nearest
rem  sea-floor low. On today's states Earth's profile stands 630-830 m above the planet's floor within 200 km of the edge;
rem  S's room is ~20 Mkm3 of elevation (wedge_room.py), up to ~55 m of sea if it filled.
rem  Arms (2 x 6 seeds = 12 runs, about 40 min); folders margin_wedge-...:
rem   False  today's PARAMS (control; must equal ocean_depth_rule-earth of sweep-20261113-201653: identity.txt)
rem   True   S on (the candidate)
rem  Rule. Vetoes on the True arm: V1 pytest passes (incl. test_margin_wedge); V2' elevation gates +-300 m on 6 of 6 seeds;
rem   V3 ledgers close (crust_budget residuals of the control's size); V4 plate_pieces_review T1 detached area <= 0.01 on
rem   every seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: True if its arm passes every veto (Earth's rise, built in front of its source to its measured profile,
rem   replaces ponding in the nearest low; a worse score becomes the next thing to work on).
rem  Cloud pre-check, seed 1014 (unused), full length (S adds sediment every step; journal 102's rule), control and S,
rem   read only, not pooled (precheck_read_1014.txt): bounded and inert on this seed. S booked 0.000 km3/yr over 250 Myr:
rem   this seed sends 0.008 km3/yr past its shelves (the six sweep seeds 0.13); at 0 Myr the sea -512 m and land p50
rem   897 / 898 m on both, highest node 5,432 / 5,433 m, gates young ocean -9 m, thin crust -41 m on both.
rem  Expected (written before the runs, not a band): the deep-sea supply is small (0.13 km3/yr over the run, 0.44 over the
rem   last 40 Myr; the shortfall is supply, not route: sediment_by_height.py), so the sea at most +55 m (S's room full),
rem   likely +10-30 m, inside the seed spread (paired sd 44-126 m of sea, 71-171 m of land p50); land p50 unresolved; the
rem   floor 0-300 km from passive edges toward Earth's profile (margin_wedge.txt); S1-S4 reported.
rem  Scored, reported (no veto): three_factor, basin, continent, hyps, column_split per arm; lith_crust2; margin_wedge
rem   (S1-S4, the floor by distance from passive edges, the booked wedge); sediment_where, sediment_by_height and
rem   wedge_room on both arms; land fraction and p25 / p50 / p90, C1, crust bands, export, volume, K1-K4, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
set P83=%CAL%\precheck-p83-sediment
%PY% -u -m scripts.run_sweep --name margin-rise --grid margin_wedge=False,True --seeds 1001,1002,1003,1004,1005,1006 --control margin_wedge-False --identity %CAL%\sweep-20261113-201653\ocean_depth_rule-earth || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P80%\lith_crust2.py "%SW%\margin_wedge-False" "%SW%\margin_wedge-True" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% -u -m planet_reviews.margin_wedge "%SW%\margin_wedge-False" "%SW%\margin_wedge-True" --json "%SW%\reviews\margin_wedge.json" > "%SW%\reviews\margin_wedge.txt" 2>&1 || goto fail
for %%a in (False True) do (
  %PY% %P82%\three_factor.py "%SW%\margin_wedge-%%a" > "%SW%\reviews\three_factor_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\basin.py "%SW%\margin_wedge-%%a" > "%SW%\reviews\basin_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\continent.py "%SW%\margin_wedge-%%a" > "%SW%\reviews\continent_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\hyps.py "%SW%\margin_wedge-%%a" > "%SW%\reviews\hyps_%%a.txt" 2>&1 || goto fail
  %PY% %P75%\column_split.py "%SW%\margin_wedge-%%a" > "%SW%\reviews\column_split_%%a.txt" 2>&1 || goto fail
)
%PY% %P83%\sediment_where.py "%SW%\margin_wedge-False" "%SW%\margin_wedge-True" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\margin_wedge-False" "%SW%\margin_wedge-True" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
%PY% %P83%\wedge_room.py "%SW%\margin_wedge-False" "%SW%\margin_wedge-True" > "%SW%\reviews\wedge_room.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
