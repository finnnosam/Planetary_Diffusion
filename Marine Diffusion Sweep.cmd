@echo off
rem Marine Diffusion Sweep (20 Nov; journal 130, ninety-seventh pass). Rule fixed 20 Nov before the runs.
rem  Why: STATE section 0, the continents' level on the datum, the spread's low side. start_cover 'accommodation' was
rem  adopted this pass (Start Cover Sweep, sweep-20261120-212152): the continents under the datum 0.245 -> 0.197 (Earth
rem  0.112), the drowned band -1,000 ... -200 m 7.9 -> 6.6% of the surface (Earth 2.0), land 26.0 -> 26.6%; the start's
rem  placement fades by 0 Myr (partial r(S, b | Tc) -0.63 at -250 Myr, -0.40 at 0 against the control's -0.37). Section
rem  0 names the next item: transport in drowned continent below wave base. The drowned continent is not starved (its
rem  cover 3,026 m against Earth's 2,611 on drowned continent); it is shaped. On the 46 km grid the landscape kernel's
rem  marine diffusion (Kd_marine 2e3 m2/yr, tuned in the 11.5 km test bed) acts on every sea node shallower than 1,000 m
rem  at every depth alike: ~400 m/Myr across a 450 m step between neighbours, below wave base as above it, and it cuts the
rem  drowned continent's crust where its cover is gone. The planet's submerged continent deepens with distance from land
rem  as a diffusive ramp, p50 -150 / -249 / -356 / -467 / -585 m at 0-50 / -100 / -150 / -200 / -300 km, where Earth's is
rem  a shelf, -51 / -147 / -212 / -236 / -293 (p75 -17 ... -89 to 300 km; precheck-p97-tail\shelf_profile.py, ETOPO and
rem  CRUST1.0 on the grid, ice out). Stratigraphic models make marine transport decay with depth (wave energy; Kaufman et
rem  al. 1991), as the plate layer's own shelf law does (shelf_kappa_m2_yr x exp(-depth / 50 m)).
rem  marine_diffusion (built this pass; tests\test_marine_diffusion.py): 'uniform' as built; 'none': the kernel's uniform
rem  term off on the plate grid; marine transport is the shelf law and the gravity paths (deltas and graded shelves, the
rem  margin wedge, fill_and_spill's turbidites). No new number. Default path identical step by step (n = 48, 70 steps).
rem  Arms (2 x 6 seeds = 12 runs), marine_diffusion:
rem   uniform  control (today's PARAMS; must equal sweep-20261120-212152's start_cover-accommodation)
rem   none     the scale-tuned uniform marine diffusion off on the plate grid
rem  Rule. Vetoes on 'none': V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close: on every
rem   seed the series' crust ledger closes (v3_split closure <= 0.001 km3/yr) and crust_budget's crust residual is the
rem   control's size by seed; a seed outside the control's range passes if v3_split puts the excess in merges / pieces /
rem   holes, or its update part, as a share of crystalline erosion, lies within this sweep's control seeds' range widened
rem   by 2 sd of them on each side; V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the
rem   highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: 'none' is adopted if it passes every veto (a construction fix: a test-bed tuned term whose depth-blind
rem   transport has no counterpart below wave base is taken off the plate grid; no number added). Fails a veto: the
rem   control kept, the failing veto is next.
rem  Design check (read, not a veto): shelf_profile's submerged-continent p50 at 0-50 and 50-100 km from land shallower
rem   than the control's (Earth -51 / -147 m), and drowned_band's continent at -1,000 ... -200 m at 0 Myr.
rem  Expected (not a band): shelves shallower near land; the drowned continent no longer cut below wave base; the risk
rem   (the pre-check): the sediment the term spread offshore stays nearer the coast or rides the gravity paths to the deep
rem   floor, so the drowned interior far from land need not rise; the sea and the young floor move with it.
rem  Scored, reported (no veto): as the Start Cover Sweep (level_split, rest_split, continent, three_factor, basin, hyps,
rem   column_split per arm; coast_profile / coast_budget, lith_crust2, thin_where, regime_split, sediment_where,
rem   sediment_by_height, thick_tail, crust_budget, v3_per_seed, v3_events, v3_split, the margins band,
rem   passive_survival, cover_accommodation, drowned_band at 0 Myr), plus shelf_profile and drowned_setting.
rem  Full-length cloud pre-check, seed 1023 (unused), control and 'none', read only, not pooled (journal 130 sec. 5;
rem   precheck-p97-tail\precheck_1023.txt): every veto reading passes: gates pooled young-ocean -111 / thin-crust -1 m
rem   (control -119 / +14); highest node in the snapshots 6,648 m (5,926); detached 0.0017 (0.0018); crust residual whole
rem   run -0.11 km3/yr (-0.14). At 0 Myr: land 0.291 (0.271); submerged continent p50 at 0-50 / 50-100 km -146 / -225 m
rem   (-187 / -272); drowned band 0.053 (0.067); but continents under the datum 0.326 (0.318), land p50 491 m (426), the
rem   sea -193 m on the datum (-110) and the ocean older (mean 95 Myr against 78). One seed; the sweep decides.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
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
set A=marine_diffusion-uniform
set B=marine_diffusion-none
%PY% -u -m scripts.run_sweep --name marine-diffusion --grid marine_diffusion=uniform,none --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261120-212152\start_cover-accommodation || goto fail
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
%PY% %P95%\v3_per_seed.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\v3_per_seed.txt" 2>&1 || goto fail
%PY% %P95%\v3_events.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\v3_events.txt" 2>&1 || goto fail
%PY% %P95%\v3_split.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\v3_split.txt" 2>&1 || goto fail
%PY% %P96%\cover_accommodation.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\cover_accommodation.txt" 2>&1 || goto fail
%PY% %P96%\drowned_band.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\drowned_band.txt" 2>&1 || goto fail
%PY% %P97%\shelf_profile.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\shelf_profile.txt" 2>&1 || goto fail
%PY% %P97%\drowned_setting.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\drowned_setting.txt" 2>&1 || goto fail
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
