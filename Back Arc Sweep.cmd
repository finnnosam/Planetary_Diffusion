@echo off
rem Back Arc Sweep (18 Nov; journal 119, ninety-second pass). Rule fixed 18 Nov before the runs.
rem  Why: STATE section 0, order item 2, the margins and the low continents. The continents' spread on the datum is
rem  134% of the land's excess, on the low side; the thin crust (total 20-30 km) is 0.53x Earth's, and pass 91
rem  (journal 118) found it short inland: 200-300 / 300-500 / 500-800 km from the edge 0.037 / 0.016 / 0.010 of the
rem  continent against Earth's 0.17 / 0.098 / 0.058. Earth's inland thin crust lies near subduction (57% within
rem  1,000 km against 19% of its band): back-arc extension, which margin_shortening's max(0, u - rollback) discards.
rem  backarc_extension (built this pass; tests\test_backarc.py): where an overriding continental trench node has
rem  u - rollback < 0 (Heuret & Lallemand 2005; Lallemand et al. 2008), the arc-forearc block steps trenchward at
rem  e = rollback - u (the plate claims the ocean node in front of the trench at e x dt x trench length per node area)
rem  and the same area opens behind the arc as pure shear of the segment's back-arc band (Earth's excess thin crust by
rem  distance to subduction, 200-1,500 km, mean 745 km: precheck-p92-backarc\thin_by_trench.txt); volume conserved
rem  (ledger 'backarc'), never below breakup_thick_m, the lithosphere stretched with the crust.
rem  The trench law decides where the branch acts (journal 119 sec. 2; precheck-p92-backarc\trench_kinematics.txt, on
rem  today's control, -40 ... 0 Myr): on the planet's continental overriding trench, beyond +-5 mm/yr, shortening /
rem  extending 0.15 / 0.62 at rollback 17, 0.32 / 0.32 at 9, 0.50 / 0.24 under Lallemand's law; Earth's continental
rem  upper plates 0.21 / 0.11 (Lallemand et al. 2008 Table 1). rollback_rule (built this pass): 'fixed' (rollback_mm_yr)
rem  or 'lallemand' (rollback = 14 - 0.24 Vsub mm/yr, their neutral-arc line in the SB04 frame); both branches.
rem  Arms (4 x 6 seeds = 24 runs), backarc_extension : rollback_rule : rollback_mm_yr
rem   False : fixed : 17      control (today's PARAMS; must equal start_sediment-all of sweep-20261117-195758)
rem   True  : fixed : 17      the extensional branch at Oncken et al. 2006's rollback (the shortening's own constant)
rem   True  : fixed : 9       the same at Earth's mean trench rollback in the SB04 frame (Lallemand et al. 2008 Table 2:
rem                           0.9 +- 3.0 cm/yr over 166 transects); it also raises the shortening (one identity)
rem   True  : lallemand : -   the trench law of Lallemand et al. 2008 (rollback falls with the subducting plate's speed)
rem  Rule. Vetoes on each back-arc arm: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close
rem   (crust_budget residuals of the control's size); V4 plate_pieces_review T1 detached area <= 0.01 on every seed;
rem   V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption by mechanism: a back-arc arm that passes every veto is adopted (it restores the measured identity's
rem   discarded half). If both pass, the one whose regime split on continental overriding trench over the last 100 Myr
rem   (regime_split.txt, six-seed means: shares extending / shortening beyond 5 mm/yr) lies nearest Earth's (Lallemand
rem   2008 Table 1, continental upper plates: extension 0.11, compression 0.21; distance = sum of the two absolute
rem   differences; ties to the earlier arm in the list).
rem   None passes: the control kept, the failing veto is the next thing to work on.
rem  Full-length cloud pre-check, seed 1019 (unused), control and the three back-arc arms, read only, not pooled
rem   (journal 119 sec. 5; precheck-p92-backarc\precheck_read_1019.txt, regime_split_1019.txt, thin_1019.txt,
rem   gate_1019.txt): every veto reading passes on all three (highest node 6,171-6,424 m against the control's 6,067;
rem   detached area <= 0.0021; gates pooled young-ocean -59 / -108 / -171, thin-crust -151 / -225 / -287 m: watch
rem   the thin-crust gate on 'lallemand'). Thin share (window) 0.24 / 0.16 / 0.21 against the control's 0.09 (Earth
rem   0.155); land fraction 0.20-0.21 against 0.26 (Earth 0.28); regime split, last 100 Myr, extending / shortening
rem   beyond 5 mm/yr 0.79 / 0.02, 0.31 / 0.26, 0.25 / 0.43 (Earth 0.11 / 0.21).
rem  Expected (not a band): thin crust 200-800 km from the edge up, beside subduction; the interior lowland up; the
rem   continents' low side lower before the sea answers; trench consumption and arc growth up a little.
rem  Scored, reported (no veto): thin_where (thin share by distance from the edge; Earth 0.17 / 0.098 / 0.058 at
rem   200-300 / 300-500 / 500-800 km), thin_setting_planet (inland thin within 1,000 km of subduction; Earth 0.57),
rem   regime_split, coast_profile / coast_budget (land below 200 m beyond 300 km inland 19.2%, rise 386 m, C1 -524 m),
rem   three_factor, basin, continent, hyps, column_split per arm, K4, the shelf, the margins band, sediment_where,
rem   thick_tail, lith_crust2, crust_budget, t3_origin; the passive record by birth type on the control
rem   (passive_survival).
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
set A=backarc_extension-False_rollback_rule-fixed_rollback_mm_yr-17
set B=backarc_extension-True_rollback_rule-fixed_rollback_mm_yr-17
set C=backarc_extension-True_rollback_rule-fixed_rollback_mm_yr-9
set D=backarc_extension-True_rollback_rule-lallemand_rollback_mm_yr-17
%PY% -u -m scripts.run_sweep --name back-arc --grid backarc_extension,rollback_rule,rollback_mm_yr=False:fixed:17,True:fixed:17,True:fixed:9,True:lallemand:17 --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261117-195758\start_sediment-all || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P92%\regime_split.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\regime_split.txt" 2>&1 || goto fail
%PY% %P92%\thin_setting_planet.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\thin_setting_planet.txt" 2>&1 || goto fail
%PY% %P91%\thin_where.py t+0000.0 "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\thin_where.txt" 2>&1 || goto fail
%PY% %P80%\lith_crust2.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
for %%a in (A B C D) do call :perarm %%a || goto fail
%PY% %P83%\sediment_where.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
%PY% %P89%\t3_origin.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\t3_origin.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:perarm
call set "ARM=%%%1%%"
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
