@echo off
rem Transfer Sweep (27 Nov; journal 153, pass 112). Rule fixed 27 Nov before the runs.
rem  Why: C (remap_rule 'overlap' + rift_area 'plates') passed every veto in the Rift Area Sweep
rem   (sweep-20261126-161947) but, made the default, failed the suite: oceanic basement thinned below 7 km with nothing
rem   booked (0.26 Mkm3 at 0 Myr, n = 217; journal 150). The review's trace and the user's 250 Myr trace (journal 151)
rem   put it in the return: a forearc cut on a world cell mixing arc and basement spread over every contributor (55% of
rem   crossings, all of the onset), and the per-plate uniform correction moving untouched basement (45%, 66% late).
rem  The repair (built this pass; composition._tent_return; tests\test_return_material.py): each world cell's weights
rem   renormalised over its live contributors; gains and losses split before mixing; a loss debited from each
rem   contributor in proportion to its material above its floor (oceanic basement keeps 7 km, the forearc law's own
rem   floor), at most all of it; gains and losses closed per plate and crust type by rescaling only the nodes that took
rem   them (factor near one); a gain with no live contributor or a loss beyond the material booked. No number added.
rem  First launch (27 Nov 17:25, sweep-20261127-173125-void): the first version closed the whole world volume and on
rem   seed 1004 at -134 Myr divided a hidden node's 102 m gain by denormal live gains (factor infinite, the run NaN);
rem   the sweep's parent process also ended unseen at ~17:54 (tab gone, no reviews). Void; the fix is the live
rem   renormalisation above (nan_probe.py reproduced and located it), and this sweep is rerun from the start.
rem  Pre-checks on the fixed version (precheck-p112-transfer, rule.md written first): unit tests (arc beside basement,
rem   ownership change, newly exposed node, the seed 1004 case) pass; trace.py seeds 1001 / 1002 10 steps and 1001 /
rem   1004 250 steps: no oceanic node below 7,000 m, no unbooked crossing, all finite; suite with C as the default:
rem   only the two default pins fail; default path identical step by step at threads 1 and 2 (71 of 71, n = 48).
rem  Arms, remap_rule and rift_area together:
rem   A nearest, world   control (today's PARAMS). Development seeds: its first seed re-run and compared with
rem                      sweep-20261122-150810\combo-01, the other five hard-linked. Confirmation seeds: run in full.
rem   C overlap, plates  the remap with the rift's opened area and the repaired return
rem  Two calls: development seeds 1001-1006 (read first), then confirmation seeds 2001-2006, whose prediction is frozen
rem   in precheck-p112-transfer\prediction.md (written before this sweep started; read only after the development seeds).
rem  Rule. Vetoes on C, on both seed sets: V1 pytest passes (the planet tests; with C as the default checked in the
rem   pass); V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close (v3_split closure <= 0.001 km3/yr; crust_budget's
rem   residual the control's size by seed, or its update share inside the control seeds' range +- 2 sd); V4 detached
rem   area <= 0.01 on every seed; V5 highest node <= 8,850 m or no higher than the control's; and the invariant the repair
rem   restores: no oceanic node below 6,999 m in any snapshot (ocean_thin.py); no RuntimeWarning or nan in any run.log.
rem  Adoption (as pass 111's rule): a numerics fix. C becomes PARAMS if it passes every veto on both seed sets; fails
rem   one: the control kept, the failing veto next. The terrain claim is made only if the frozen prediction holds on
rem   2001-2006 (all seven deciding lines); a failed line is reported and is the next problem.
rem  Scored, reported: geography (both seed sets, window and endpoint); the continents on the datum (continent,
rem   level_split, three_factor), hypsometry; run time per arm.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set NUMBA_NUM_THREADS=2
set P82=%CAL%\precheck-p82-stepback
set P94=%CAL%\precheck-p94-level
set P95=%CAL%\precheck-p95-v3
set P112=%CAL%\precheck-p112-transfer
set A=remap_rule-nearest_rift_area-world
set C=remap_rule-overlap_rift_area-plates
%PY% -u -m scripts.run_sweep --name transfer --grid remap_rule,rift_area=nearest:world,overlap:plates --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261122-150810\combo-01 --borrow-control || goto fail
call :latest
call :probes || goto fail
set "SWDEV=%SW%"
set PLANET_SEEDS=2001,2002,2003,2004,2005,2006
%PY% -u -m scripts.run_sweep --name transfer-confirm --tests tests\test_return_material.py --grid remap_rule,rift_area=nearest:world,overlap:plates --seeds 2001,2002,2003,2004,2005,2006 --control %A% || goto fail
call :latest
call :probes || goto fail
echo Done: development %SWDEV%, confirmation %SW%
pause
exit /b 0
:latest
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & exit /b 0)
exit /b 0
:probes
echo Probes on %SW%
%PY% %P112%\ocean_thin.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\ocean_thin.txt" 2>&1 || exit /b 1
%PY% %P95%\v3_split.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\v3_split.txt" 2>&1 || exit /b 1
%PY% %P95%\v3_per_seed.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\v3_per_seed.txt" 2>&1 || exit /b 1
for %%a in (A C) do call :perarm %%a || exit /b 1
exit /b 0
:perarm
call set "ARM=%%%1%%"
%PY% %P94%\level_split.py "%SW%\%ARM%" > "%SW%\reviews\level_split_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\three_factor.py "%SW%\%ARM%" > "%SW%\reviews\three_factor_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\continent.py "%SW%\%ARM%" > "%SW%\reviews\continent_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\hyps.py "%SW%\%ARM%" > "%SW%\reviews\hyps_%1.txt" 2>&1 || exit /b 1
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
