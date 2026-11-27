@echo off
rem Pelagic Probes (8 Nov; journal 99, seventy-fifth pass). Finishes the Pelagic Sweep (sweep-20261108-150308): its
rem  plate_pieces_review stopped on the 1.1 arm (seed 1002, -60 Myr: one node still labelled with plate 17, dropped that
rem  step when its raster emptied; trench_origin now skips such a label), so run_sweep failed and the .cmd stopped
rem  before its probes. This reruns the plate review (veto V4: T1 detached area <= 0.01 on every seed of the 1.1 arm)
rem  and the probes, each into the sweep's reviews\. No runs; a few minutes. The rule is the Pelagic Sweep's, unchanged.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P67=%CAL%\precheck-p67-survival
set P69=%CAL%\precheck-p69-gate
set P71=%CAL%\precheck-p71-window
set P74=%CAL%\precheck-p74-pelagic
set SW=%CAL%\sweep-20261108-150308
set A=%SW%\pelagic_m_myr-0
set B=%SW%\pelagic_m_myr-1.1
echo Plate review on %SW%
%PY% -u -m scripts.plate_pieces_review "%A%" "%B%" --control pelagic_m_myr-0 > "%SW%\reviews\plate_pieces_review.txt" 2>&1 || goto fail
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
echo FAILED, see the file just written in %SW%\reviews
pause
exit /b 1
