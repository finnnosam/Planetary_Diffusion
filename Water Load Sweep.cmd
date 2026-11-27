@echo off
rem Water Load Sweep (19 Nov; journal 124, ninety-fourth pass). Rule fixed 19 Nov before the runs.
rem  Why: STATE section 0, the continents' level on the datum. Pass 94 split the level by crust bin and term, planet
rem  (today's control) against Earth on the same nodes, alike statistics (precheck-p94-level\level_split.txt): over
rem  30-45 km crust the planet stands 272 m low; the largest term is not a candidate but the column's remainder
rem  (rest -132 m; -627 m on 30-35 km). Its split (rest_split.txt): surface.elevation water-loads a thinned crust's own
rem  Airy depth wherever it lies below sea level (crust x rw / (rm - rw)), also where the column's sediment fills it
rem  above sea level and is booked as a dry load: no water is in that column. -306 m on 30-35 km crust (72% of the
rem  bin), -200 m on the continent's mean. On a law with one water decision per column (1-D Airy; the law
rem  elevation_check's thin-crust gate already reads, which is why that gate reads -173 ... -230 m on today's model)
rem  the planet's 30-35 km crust stands 327 m below its own column, Earth's 314 m above it.
rem  water_load_rule (built this pass; tests\test_water_load.py): 'parts' as built; 'column' decides once per
rem  continental node by its air-loaded local level (crust, densified root, sediment's bulk height less its dry mass)
rem  against sea level: wet columns water-load crust and sediment, dry ones neither. A correctness fix, no value.
rem  Default path identical step by step (n = 48, 70 steps, threads 1).
rem  Arms (2 x 6 seeds = 12 runs), water_load_rule:
rem   parts   control (today's PARAMS; must equal sweep-20261119-182758's backarc_extension-True_..._mm_yr-9)
rem   column  one water decision per column
rem  Rule. Vetoes on 'column': V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close
rem   (crust_budget residuals of the control's size); V4 plate_pieces_review T1 detached area <= 0.01 on every seed;
rem   V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: 'column' is adopted if it passes every veto (a correctness fix: the model's law made consistent with
rem   itself and with its own gate). Fails a veto: the control kept, the failing veto is the next thing to work on.
rem  Design check (read, not a veto): the pooled thin-crust gate offset nearer zero than the control's on >= 5 of 6
rem   seeds (the gate reads the column law 'column' builds; flexure stays in it).
rem  Expected (not a band): continents on the datum up on 20-40 km crust; drowned continental crust and the
rem   -1,000 ... -200 m band down; land fraction up; the sea up a little (columns that rise hold less water).
rem  Scored, reported (no veto): level_split and rest_split per arm (precheck-p94-level), continent, three_factor,
rem   basin, hyps, column_split per arm, coast_profile / coast_budget, lith_crust2, thin_where, regime_split,
rem   sediment_where, sediment_by_height, thick_tail, crust_budget, the margins band, passive_survival.
rem  Full-length cloud pre-check, seed 1007 (unused), control and 'column', read only, not pooled (journal 124
rem   sec. 6): every veto reading passes on 'column': highest node 6,418 m (control 6,421); detached area 0.0034
rem   (0.0014); gates pooled young-ocean -114 / thin-crust +5 m (control -51 / -145: the design check holds);
rem   whole-run ledger residuals +0.11 / -0.12 km3/yr (control +0.09 / -0.10), last 40 Myr continental +0.72
rem   (control +0.23; the six-seed control -0.09): read V3 per seed. The runs part early (48 plates at 0 Myr
rem   against 28). Continents on the datum p25 / p50 / p75 -283 / 106 / 561 m, under the datum 0.37 (control
rem   -393 / 48 / 465, 0.45); the sea +39 m on the datum (control -38): land 0.192 (0.207), window 0.217 (0.235).
rem   The continents rise and the sea rises with them; one seed cannot split the land fraction from the drift.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
set P83=%CAL%\precheck-p83-sediment
set P84=%CAL%\precheck-p84-thicktail
set P89=%CAL%\precheck-p89-startsed
set P91=%CAL%\precheck-p91-thin
set P92=%CAL%\precheck-p92-backarc
set P94=%CAL%\precheck-p94-level
set A=water_load_rule-parts
set B=water_load_rule-column
%PY% -u -m scripts.run_sweep --name water-load --grid water_load_rule=parts,column --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261119-182758\backarc_extension-True_rollback_rule-fixed_rollback_mm_yr-9 || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P92%\regime_split.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\regime_split.txt" 2>&1 || goto fail
%PY% %P91%\thin_where.py t+0000.0 "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\thin_where.txt" 2>&1 || goto fail
%PY% %P80%\lith_crust2.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
for %%a in (A B) do call :perarm %%a || goto fail
%PY% %P83%\sediment_where.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:perarm
call set "ARM=%%%1%%"
%PY% %P94%\level_split.py "%SW%\%ARM%" > "%SW%\reviews\level_split_%1.txt" 2>&1 || exit /b 1
%PY% %P94%\rest_split.py "%SW%\%ARM%" > "%SW%\reviews\rest_split_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\three_factor.py "%SW%\%ARM%" > "%SW%\reviews\three_factor_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\basin.py "%SW%\%ARM%" > "%SW%\reviews\basin_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\continent.py "%SW%\%ARM%" > "%SW%\reviews\continent_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\hyps.py "%SW%\%ARM%" > "%SW%\reviews\hyps_%1.txt" 2>&1 || exit /b 1
%PY% %P75%\column_split.py "%SW%\%ARM%" > "%SW%\reviews\column_split_%1.txt" 2>&1 || exit /b 1
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
