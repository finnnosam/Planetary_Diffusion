@echo off
rem Rift Area Sweep (26 Nov; journal 148, pass 110). Rule fixed 26 Nov before the runs.
rem  Why: the Remap Sweep (sweep-20261126-141105) failed V3 on remap_rule 'overlap': the crust budget's update part,
rem   whole run, -0.36 km3/yr against the control's -0.08 (update / crystalline erosion -0.123 ... -0.177 on six of six
rem   seeds against the control's -0.020 ... -0.075; +-2 sd -0.083 ... +0.003); V1, V2', V4, V5 passed. Its geography
rem   moved toward Earth on 1001-1006 (window: holes 358 -> 221, low flat interior 3.8 -> 6.8% of land, low land among
rem   far land 3.6 -> 10.3%, crossings 53.4 -> 44.5; Earth 13 / 14.2 / 19.4 / 19.1); the continents sat lower on the
rem   datum (share under it 0.193 -> 0.273, Earth 0.112).
rem  The leak (precheck-p110-rift, rule-rift.md committed first, md5 35594f2d; cloud seed 1028, warm state at -20 Myr,
rem   10 steps per arm): the world's pure shear strains the zone by each gap's world cell; the plates switch on the
rem   nodes whose forward nearest node is the gap, on average 0.91 of that area (13% of gaps open nothing, ~half
rem   more than their cell). The new plate nodes take 0.47 km3/yr less crust than the zone gives, the same under
rem   'nearest' (-0.472) and 'overlap' (-0.468): reading R1. 'nearest''s writeback skipped ~7% of the thinning, which
rem   hid it.
rem  rift_area (built this pass; tests\test_rift_area.py): 'world' (today) / 'plates': the strain is the area the
rem   owning plate's raster switches on over each gap this step, so the zone gives exactly the volume the new plate
rem   nodes take (probe: thinned = placed to rounding); a gap no plate opens over takes no volume and is shown as the
rem   zone's crust. No number added. Default path identical step by step at threads 1 and 2 (70 of 70 states each,
rem   n = 48), a run against itself at threads 2 identical.
rem  Full-length cloud pre-check, seed 1032 (unused), control and C, threads 1, read only, not pooled
rem   (precheck-p110-rift\precheck_1032): no numerical warnings; gates -127 / +28 m (control -110 / +30), highest
rem   node 6,749 m (6,051), detached 0.0021 (0.0022); V3 update / crystalline erosion -0.043 (control +0.021; overlap
rem   without the opened area -0.146 on seed 1031), residual +0.23 (+0.16); crust -1.7% (the control -7.5% on this
rem   seed). Window -40 ... 0: holes 289 -> 221, low flat interior 2.9 -> 5.1%, low land among far land 1.9 -> 7.6%,
rem   crossings 41.1 -> 44.6, land 26.7 -> 26.7%. Continents on the datum lower: share under it 0.198 -> 0.361.
rem  Arms (borrowed control: 1 + 6 runs), remap_rule and rift_area together:
rem   A nearest, world   control (today's PARAMS); its first seed re-run and compared with
rem                      sweep-20261122-150810\combo-01, the other five hard-linked from it
rem   C overlap, plates  the remap with the rift's opened area
rem   (the Remap Sweep's arm overlap, world is read beside them for attribution: same code, deterministic)
rem  Rule. Vetoes on C: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close (the series'
rem   crust ledger closes, v3_split closure <= 0.001 km3/yr; crust_budget's residual the control's size by seed, or
rem   v3_split puts the excess in merges / pieces / holes or its update share within this sweep's control seeds' range
rem   +- 2 sd); V4 T1 detached area <= 0.01 on every seed; V5 on every seed the highest node over the run <= 8,850 m or
rem   no higher than the control's.
rem  Adoption: a numerics fix (the plate rasters carry what the world shows, and the rift's new crust is the volume the
rem   zone gives). C is adopted (remap_rule 'overlap' and rift_area 'plates' become PARAMS) if it passes every veto.
rem   Fails one: the control kept, the failing veto next. The terrain claim waits for the confirmation seeds 2001-2006
rem   against a prediction frozen in the journal before that .cmd is handed over.
rem  Scored, reported: geography (STATE 4a) on all three, endpoint and window; the continents on the datum (continent,
rem   three_factor, hyps) first among the rest; the Rift Traction Sweep's reviews and probes; run time per arm.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set NUMBA_NUM_THREADS=2
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
set P83=%CAL%\precheck-p83-sediment
set P84=%CAL%\precheck-p84-thicktail
set P91=%CAL%\precheck-p91-thin
set P92=%CAL%\precheck-p92-backarc
set P94=%CAL%\precheck-p94-level
set P95=%CAL%\precheck-p95-v3
set P96=%CAL%\precheck-p96-spread
set P97=%CAL%\precheck-p97-tail
set P102=%CAL%\precheck-p102-share
set P103=%CAL%\precheck-p103-meet
set P107=%CAL%\precheck-p107-transfer
set A=remap_rule-nearest_rift_area-world
set C=remap_rule-overlap_rift_area-plates
set D=%CAL%\sweep-20261126-141105\remap_rule-overlap
%PY% -u -m scripts.run_sweep --name riftarea --grid remap_rule,rift_area=nearest:world,overlap:plates --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261122-150810\combo-01 --borrow-control || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P107%\geography.py "%SW%\%A%" "%SW%\%C%" "%D%" --json "%SW%\reviews\geography.json" > "%SW%\reviews\geography.txt" 2>&1 || goto fail
%PY% %P95%\v3_split.py "%SW%\%A%" "%SW%\%C%" "%D%" > "%SW%\reviews\v3_split.txt" 2>&1 || goto fail
%PY% %P95%\v3_per_seed.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\v3_per_seed.txt" 2>&1 || goto fail
%PY% %P95%\v3_events.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\v3_events.txt" 2>&1 || goto fail
for %%a in (A C) do call :perarm %%a || goto fail
%PY% -u -m planet_reviews.dynamic_topography "%SW%\%A%" "%SW%\%C%" --json "%SW%\reviews\dynamic_topography.json" > "%SW%\reviews\dynamic_topography.txt" 2>&1 || goto fail
%PY% %P92%\regime_split.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\regime_split.txt" 2>&1 || goto fail
%PY% %P91%\thin_where.py t+0000.0 "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\thin_where.txt" 2>&1 || goto fail
%PY% %P80%\lith_crust2.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
%PY% %P83%\sediment_where.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
%PY% %P96%\drowned_band.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\drowned_band.txt" 2>&1 || goto fail
%PY% %P96%\dyn_tail.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\dyn_tail.txt" 2>&1 || goto fail
%PY% %P97%\shelf_profile.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\shelf_profile.txt" 2>&1 || goto fail
%PY% %P102%\plate_continent.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\plate_continent.txt" 2>&1 || goto fail
%PY% %P103%\cont_spectrum.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\cont_spectrum.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:perarm
call set "ARM=%%%1%%"
%PY% %P94%\level_split.py "%SW%\%ARM%" > "%SW%\reviews\level_split_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\three_factor.py "%SW%\%ARM%" > "%SW%\reviews\three_factor_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\basin.py "%SW%\%ARM%" > "%SW%\reviews\basin_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\continent.py "%SW%\%ARM%" > "%SW%\reviews\continent_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\hyps.py "%SW%\%ARM%" > "%SW%\reviews\hyps_%1.txt" 2>&1 || exit /b 1
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
