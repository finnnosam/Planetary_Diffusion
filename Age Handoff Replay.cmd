@echo off
rem Age Handoff Replay (23 Nov; journal 142, pass 105). The review (journal 141 section 1 item 4) found that the
rem  landscape stage gave active belts old-rock erodibility (0.15 instead of 1.0) for its whole 30 Myr window. Pass 105
rem  fixed the hand-off (coupling.erodibility_map, rule 'plate': the plate layer's own law on its nodes, then sampled).
rem  This replays the review's Chain run (seed 1001, every-Myr snapshots) through the landscape and the meso twice on
rem  the Chain run's own forcing (read, not rewritten): once with the old rule ('legacy', which must reproduce the Chain
rem  run's landscape exactly) and once with the fix ('plate'). The plate layer is untouched (no plate code changed).
rem  Measures fixed in journal 142 before the runs; read in pass 106. Not scored, not a sweep. About 45 minutes.
rem  Output: outputs\planet-landscape\age-replay-20261123 (legacy\, plate\, compare.txt).
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set NUMBA_NUM_THREADS=2
set CH=outputs\review-astra\runs\Chain-360f894f
set RP=outputs\planet-landscape\age-replay-20261123
%PY% -m pytest -q tests\test_coupling.py || goto fail
%PY% -u -m scripts.run_planet_landscape --plates "%~dp0%CH%\plates\seed1001" --forcing %CH%\landscape\forcing --out %RP%\legacy\landscape --age-rule legacy || goto fail
%PY% -u -m scripts.planet_meso_check --landscape %RP%\legacy\landscape --out %RP%\legacy\meso || goto fail
%PY% -u -m scripts.run_planet_landscape --plates "%~dp0%CH%\plates\seed1001" --forcing %CH%\landscape\forcing --out %RP%\plate\landscape --age-rule plate || goto fail
%PY% -u -m scripts.planet_meso_check --landscape %RP%\plate\landscape --out %RP%\plate\meso || goto fail
%PY% outputs\planet-calibration\precheck-p105-age\replay_compare.py %CH% %RP% > %RP%\compare.txt || goto fail
type %RP%\compare.txt
echo Done: %RP%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
