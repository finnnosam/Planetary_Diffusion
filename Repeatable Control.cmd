@echo off
rem Repeatable Control (22 Nov; journal 140, pass 104). No arm is judged: today's model (PARAMS) on six seeds, the
rem  reference the next sweep borrows its control from, and a test that a run on this machine repeats itself.
rem  Why: the Rift Traction Sweep's borrowed control did not borrow (check seed 1001 first differs at -243 Myr, 0.45 km3
rem   of continent, step 7), as the Placed Share Sweep's did not (-242 Myr). Pass 104 found the cause in the device
rem   shell: at threads 2 HEAD differs from HEAD itself at step 7 (7 of 41 states identical, n = 48), first in the
rem   mantle flow's coefficients: ducc0's adjoint transform (points to coefficients) spreads the points in an order set
rem   by the threads' timing (4 of 44 calls on the slab parcels did not repeat). The code had not changed the default
rem   path. Fixed this pass: the adjoint runs on one thread (mantle_flow.point_alm, SphereFlexure.deflection); the
rem   synthesis keeps its threads. Device shell: threads 2 against itself 71 of 71 (seed 1002) and 41 of 41 (seed 1001);
rem   threads 1 against HEAD 71 of 71 (the arithmetic at threads 1 is unchanged).
rem  Today's model: rift_traction True (adopted in pass 104 from sweep-20261121-220857). The fix changes the arithmetic
rem   at threads 2, so this run is a new trajectory of the adopted model: its identity review against
rem   sweep-20261121-220857\rift_traction-True is expected DIFFERENT and is recorded, not a veto.
rem  The test: seed 1001 is run a second time into <sweep>\repeat and compared with the sweep's own seed 1001
rem   (reviews\repeat.txt). IDENTICAL: runs here repeat; the next sweep borrows from this sweep's rift_traction-True.
rem   DIFFERENT: something else on this machine is not repeatable; the next sweep runs its control in full, and pass
rem   105 looks for it.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set NUMBA_NUM_THREADS=2
%PY% -u -m scripts.run_sweep --name repeatable-control --grid rift_traction=True --seeds 1001,1002,1003,1004,1005,1006 --control rift_traction-True --identity %CAL%\sweep-20261121-220857\rift_traction-True || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Repeat of seed 1001 into %SW%\repeat (about 20 minutes)
%PY% -u -m scripts.simulate_plates --n 217 --seeds 1001 --years 250 --out "%SW%\repeat" --param threads=2 --save-every 5 --param rift_traction=True > "%SW%\reviews\repeat_run.log" 2>&1 || goto fail
%PY% tools\compare_runs.py "%SW%\rift_traction-True\seed1001" "%SW%\repeat\seed1001" > "%SW%\reviews\repeat.txt" 2>&1
type "%SW%\reviews\repeat.txt" | findstr /b "IDENTICAL DIFFERENT"
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
