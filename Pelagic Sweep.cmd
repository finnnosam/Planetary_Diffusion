@echo off
rem Pelagic Sweep (8 Nov; journal 97, seventy-fourth pass). Rule fixed 8 Nov before the runs.
rem  The choosing reading (pass 74's ten-pass step-back, journal 97): the land is the scorecard's largest failure and
rem  rose this pass (p50 1,084 -> ~1,300 m with rift_line 'force'; Earth 384); a third of the excess is the low sea
rem  (-564 m on the datum), whose cause is a deep floor (4,574 against 4,205 m) that carries half Earth's sediment
rem  (94 against 191 Mkm3 bulk on the label-free ocean). Far from land and trenches the planet's floor is bare: 11 m at
rem  54 Myr where Earth's holds 150 m at 55 Myr (>= 1,200 km from the continent, >= 500 km from Bird's trenches): only
rem  pelagic rain (biogenic and wind-blown) reaches there, and the model has none (pelagic_m_myr 0).
rem  The value is Earth's (precheck-p74-pelagic\earth_pelagic.py, GlobSed v3 on Seton 2020 ages, n = 217): grain
rem  thickness by the ocean's Athy column, sum(a S) / sum(a age) = 1.41 / 1.12 / 1.01 m/Myr at >= 800 / 1,200 /
rem  1,600 km from the continent: 1.1 m/Myr of grain (2.7 bulk) on every oceanic node. Uniform (the carbonate
rem  compensation depth's control left out, assumed). No new fitted number.
rem  Expected (written before the runs, not a band): ~50 Mkm3 more floor sediment (bulk), the floor ~100 m shallower,
rem  the sea ~100 m higher, land p50 and land fraction down by about as much; open ocean (>= 1,200 / 500 km) toward
rem  Earth's 150 m. Passive edge and plates hardly change (a load, not a force).
rem  Pre-check (cloud seed 1007, read only, not pooled; precheck-p74-pelagic): sea -575 -> -531 m, land 32.5 -> 29.8%,
rem  ocean sediment 56 -> 137 Mkm3 bulk (Earth 191); the far ocean overshoots (274 m against 150): the floor keeps its
rem  sediment ~3x longer than it lives (journal 97 section 6). Rule unchanged.
rem  Arms (2 x 6 seeds = 12 runs, about 40-55 min); folders pelagic_m_myr-...:
rem   0    today's PARAMS (control; must equal rift_line-force of sweep-20261108-124109: identity.txt)
rem   1.1  + Earth's pelagic rain
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers
rem   close (tests, incl. test_pelagic); V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed
rem   the highest node over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: pelagic_m_myr 1.1 if its arm passes every veto (a process Earth has and the model lacks, at Earth's
rem   measured rate); not a score.
rem  Scored, reported (no veto): floor sediment stock and its ledger (crust_budget, now with the pelagic term), the
rem   open ocean by class (planet_open, ocean_sediment_probe), deep floor and sea on the datum, land fraction and land
rem   p25 / p50 / p90, coast_budget, K1-K4, edge measures, margins band, export, volume, label-free area, T2-T5, plates.
rem  Runs with scripts.run_sweep: tests, the 12 runs and the reviews under one progress line; then the probes below,
rem  each into <sweep>\reviews\*.txt.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P67=%CAL%\precheck-p67-survival
set P69=%CAL%\precheck-p69-gate
set P71=%CAL%\precheck-p71-window
set P74=%CAL%\precheck-p74-pelagic
%PY% -u -m scripts.run_sweep --name pelagic --grid pelagic_m_myr=0,1.1 --seeds 1001,1002,1003,1004,1005,1006 --control pelagic_m_myr-0 --identity %CAL%\sweep-20261108-124109\rift_line-force || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
set A=%SW%\pelagic_m_myr-0
set B=%SW%\pelagic_m_myr-1.1
echo Probes on %SW%
%PY% %P74%\planet_open.py "%A%" "%B%" > "%SW%\reviews\planet_open.txt" 2>&1 || goto fail
%PY% -u -m scripts.ocean_sediment_probe "%A%" --no-arrivals > "%SW%\reviews\ocean_sediment_probe_0.txt" 2>&1 || goto fail
%PY% -u -m scripts.ocean_sediment_probe "%B%" --no-arrivals > "%SW%\reviews\ocean_sediment_probe_1.1.txt" 2>&1 || goto fail
%PY% %P69%\edge_split.py "%A%" "%B%" > "%SW%\reviews\edge_split.txt" 2>&1 || goto fail
%PY% %P71%\edge_window.py "%A%" "%B%" > "%SW%\reviews\edge_window.txt" 2>&1 || goto fail
%PY% %P67%\margin_age_noearth.py "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\margin_age_1.1.txt" 2>&1 || goto fail
%PY% %P67%\w_read.py "%A%\seed1001" "%A%\seed1002" "%A%\seed1003" "%A%\seed1004" "%A%\seed1005" "%A%\seed1006" "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\w_read.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
