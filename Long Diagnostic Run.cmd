@echo off
rem Long Diagnostic Run (22 Nov; journal 140, pass 104). Asked by the user on 21 Nov (STATE section 5, item 8): now and
rem  then, today's model on one seed for 500 Myr, read for drift (crust volume, sediment stock, plate count, the
rem  continents' level, the slab record), not scored against Earth. Beyond 250 Myr split_rule 'earth' holds Earth's last
rem  rate (a known limit). It can run beside Repeatable Control.cmd (that one uses 12 of the 16 threads, this one 2).
rem  About an hour. Output: outputs\planet-plates\long-500-seed1001-20261122 (stats every Myr, snapshots every 10 Myr).
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set NUMBA_NUM_THREADS=2
%PY% -u -m scripts.simulate_plates --n 217 --seeds 1001 --years 500 --out outputs\planet-plates\long-500-seed1001-20261122 --param threads=2 --save-every 10 || goto fail
echo Done: outputs\planet-plates\long-500-seed1001-20261122
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
