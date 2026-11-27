@echo off
rem Remap Sweep (26 Nov; journal 145, pass 108; re-run after pass 109's fix, journal 147). Rule fixed 26 Nov before
rem  the first run; unchanged.
rem  Why: STATE section 0 item 1, the excess shoreline and the missing inland lowland. The floor's second trace (cloud
rem   seed 1028, 20 steps; outputs\planet-calibration\precheck-p108-floor, rule-floor.md committed first, md5 842142d2):
rem   - smoothing crust + sediment + thermal age at the end of a step removes 57 / 46 / 61% of the excess holes /
rem     crossings / low-flat-interior deficit; the crust alone makes the surface rougher (the sediment fills compensate
rem     node-scale crust roughness): no single field, the materials' node-scale roughness jointly (outcome O2).
rem   - the composition's jump is the plate rasters' departure from the world: the world's own end state carried to the
rem     next rotation by the same nearest-node rule adds +13 holes / -1.4 crossings; the plate rasters read there add
rem     +172 / +9.5 (validated against the next step's composition within 2%). Journal 144's suspect (a) refuted.
rem  remap_rule (built this pass; tests\test_remap_tent.py): 'nearest' (today) / 'tent' / 'overlap'. 'overlap':
rem   first-order conservative remapping; the composition takes crust, sediment, rock, thermal and formation age as the
rem   area mean of the plate cells each world cell overlaps (sampled at nine points), and a step's increments go back
rem   by the same overlaps, each plate's gains and losses per crust type conserved.
rem   No number added. Default path identical step by step at threads 1 and 2 (71 of 71 states each, n = 48).
rem  The arm's trace (same warm state; default in brackets): end of step holes 279 (416), crossings 45.8 (53.7), low
rem   flat interior 9.0% of land (5.1), low land among far land 16.1% (8.9; Earth 19.4). 'tent' (built, off, not swept)
rem   roughens its own plate rasters step by step (its weights do not sum to one per plate node); 'overlap', the
rem   conservative build of the same arm, does not.
rem  First run (sweep-20261126-093644) stopped at scoring: overlap seed 1005 went non-finite at -150 Myr. Cause (journal
rem   147): the return rescaled each plate's gains and losses separately; the remap averages erosion's node-scale cut
rem   next to a fill toward zero, and the rescale multiplied what was left by up to 1e14. Fixed: each plate's net
rem   volume per crust type conserved by an additive correction (world nodes typed by the crust they read), plate crust
rem   clamped at zero and the clamp booked (return_thick_clamped_m3). Tests 26 of 26 (new: a cut next to a fill is not
rem   amplified; crust never below zero); default path identical at threads 1 and 2 (71 of 71 each). The first run's
rem   control check seed matched sweep-20261122-150810\combo-01 on 780 of 780 fields.
rem  Full-length cloud pre-check of the fix, seed 1031 (unused), control and C, read only, not pooled
rem   (precheck-p108-floor\precheck_1031): no numerical warnings (every overlap run of the first build had them);
rem   gates -118 / +40 m (control -97 / +0), highest node 6,138 m (6,324), detached 0.0030 (0.0022), clamped crust
rem   0.0014 km3 over the run. V3 still at risk: update / crystalline erosion -0.146 (control +0.011); crust volume
rem   -5.4% (-3.3%). Suspected cause unchanged (journal 147): the exact return no longer hides the rift's stretching
rem   imbalance. Window -40 ... 0: holes 353 -> 227, crossings 49.8 -> 43.5, low flat interior 2.7 -> 3.2%, low land
rem   among far land 2.8 -> 5.7%, land 27.3 -> 25.0%. The fix on the trace: holes 257, crossings 44.1, low flat 9.6%
rem   (first build 284 / 46.0 / 8.9; control 416 / 53.7 / 5.1).
rem  Arms (borrowed control: 1 + 6 runs), remap_rule:
rem   A nearest  control (today's PARAMS); its first seed re-run and compared with sweep-20261122-150810\combo-01, the
rem              other five hard-linked from it
rem   C overlap
rem  Rule. Vetoes on C: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close (the series'
rem   crust ledger closes, v3_split closure <= 0.001 km3/yr; crust_budget's residual the control's size by seed, or
rem   v3_split puts the excess in merges / pieces / holes or its update share within this sweep's control seeds' range
rem   +- 2 sd); V4 T1 detached area <= 0.01 on every seed; V5 on every seed the highest node over the run <= 8,850 m or
rem   no higher than the control's.
rem  Adoption: the remap is a numerics fix (the plate rasters must carry what the world shows; the trace measured their
rem   departure). C is adopted if it passes every veto. Fails one: the control kept, the failing veto next. The terrain
rem   claim waits for the confirmation seeds 2001-2006 against the prediction frozen in journal 145 sec. 6 (read here on
rem   1001-1006 too, not as the confirmation).
rem  Scored, reported: geography (STATE 4a) on both arms, endpoint and window; the Rift Traction Sweep's reviews and
rem   probes; run time per arm (the remap's cost; in the cloud none).
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set NUMBA_NUM_THREADS=2
set P75=%CAL%\precheck-p75-stepback
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
set P98=%CAL%\precheck-p98-low
set P101=%CAL%\precheck-p101-degree
set P102=%CAL%\precheck-p102-share
set P103=%CAL%\precheck-p103-meet
set P107=%CAL%\precheck-p107-transfer
set A=remap_rule-nearest
set C=remap_rule-overlap
%PY% -u -m scripts.run_sweep --name remap --grid remap_rule=nearest,overlap --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261122-150810\combo-01 --borrow-control || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P107%\geography.py "%SW%\%A%" "%SW%\%C%" --json "%SW%\reviews\geography.json" > "%SW%\reviews\geography.txt" 2>&1 || goto fail
%PY% -u -m planet_reviews.dynamic_topography "%SW%\%A%" "%SW%\%C%" --json "%SW%\reviews\dynamic_topography.json" > "%SW%\reviews\dynamic_topography.txt" 2>&1 || goto fail
%PY% %P92%\regime_split.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\regime_split.txt" 2>&1 || goto fail
%PY% %P91%\thin_where.py t+0000.0 "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\thin_where.txt" 2>&1 || goto fail
%PY% %P80%\lith_crust2.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
for %%a in (A C) do call :perarm %%a || goto fail
%PY% %P83%\sediment_where.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
%PY% %P95%\v3_per_seed.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\v3_per_seed.txt" 2>&1 || goto fail
%PY% %P95%\v3_events.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\v3_events.txt" 2>&1 || goto fail
%PY% %P95%\v3_split.py "%SW%\%A%" "%SW%\%C%" > "%SW%\reviews\v3_split.txt" 2>&1 || goto fail
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
