@echo off
rem Start Cover Sweep (20 Nov; journal 129, ninety-sixth pass). Rule fixed 20 Nov before the runs.
rem  Why: STATE section 0, the continents' level on the datum. sediment_rock 'burial' was adopted this pass (Rock Cover
rem  Sweep, sweep-20261120-190819): the continents on the datum p50 349 m (Earth 328), p75 796 (744), and
rem  three_factor's level 2 m; what is left is the spread (+120 m), and its larger side is the low one: the continents'
rem  p10 / p25 -382 / -17 m (Earth -12 / 109), 24.5% under the datum (11%). The drowned continent is not Earth's kind:
rem  on the planet normal crust (median 31 km, half 30-40 km) ~400 km inland on old lithosphere, 7.9% of the surface
rem  at -1,000 ... -200 m about the sea; on Earth thin rifted margin (24 km, 82% within 300 km of the edge), 2.0%
rem  (precheck-p96-spread\drowned_band.py). It comes in with the start: at -250 Myr 5.2% of the surface is continent at
rem  -1,000 ... -200 m, crust 37 km, 3,000 km inside the supercontinent, with 296 m of cover. The start's cover is
rem  placed by a copula of the crust and a random field, blind to the surface: its partial rank correlation with the
rem  column's bare level given the crust is +0.20 at -250 Myr (cover under the datum 623 m, above 1,594), Earth's
rem  stable domain -0.76 (under 2,059, above 1,430; precheck-p96-spread\cover_accommodation.py). On Earth sediment lies
rem  where its column has room for it.
rem  start_cover (built this pass; tests\test_start_cover.py): 'copula' as built; 'accommodation': the copula of the
rem  crust and the start's bare level (Airy line of the crust before the rim is thinned, C's term of its lithosphere
rem  age, the air-loaded dynamic topography at the start), with Earth's correlations of the stable cover with both
rem  (EARTH_STABLE_SED_BARE -0.333 / -0.787, CRUST1.0 and ETOPO on the stable domain, ice sheets out) and the start's
rem  own correlation between them; the rest the same random field. Distribution, volume and every other draw of the
rem  start unchanged; no fitted number. Default path identical step by step (n = 48, 70 steps, threads 1).
rem  Arms (2 x 6 seeds = 12 runs), start_cover:
rem   copula         control (today's PARAMS; must equal sweep-20261120-190819's sediment_rock-burial)
rem   accommodation  the start's cover where its column has room for it
rem  Rule. Vetoes on 'accommodation': V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close:
rem   on every seed the series' crust ledger closes (v3_split closure <= 0.001 km3/yr) and crust_budget's crust
rem   residual is the control's size by seed; a seed outside the control's range passes if v3_split puts the excess in
rem   merges / pieces / holes, or its update part, as a share of crystalline erosion, lies within this sweep's control
rem   seeds' range widened by 2 sd of them on each side; V4 plate_pieces_review T1 detached area <= 0.01 on every seed;
rem   V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: 'accommodation' is adopted if it passes every veto (a construction fix with measured correlations: Earth's
rem   cover lies in the room its column makes). Fails a veto: the control kept, the failing veto is next.
rem  Design check (read, not a veto): cover_accommodation's partial r(S, b | Tc) at -250 Myr <= -0.5 (control about
rem   +0.2, Earth -0.76), and drowned_band's continent at -1,000 ... -200 m at 0 Myr below the control's.
rem  Expected (not a band): continents under the datum down at the start and at 0 Myr; land up; the risk (the
rem   pre-check): land p50 up (the filled lows emerge, the bare highs stay high), and the start's effect fades as the
rem   field and the plates move.
rem  Scored, reported (no veto): as the Rock Cover Sweep (level_split, rest_split, continent, three_factor, basin, hyps,
rem   column_split per arm; coast_profile / coast_budget, lith_crust2, thin_where, regime_split, sediment_where,
rem   sediment_by_height, thick_tail, crust_budget, v3_per_seed, v3_events, v3_split, the margins band,
rem   passive_survival), plus cover_accommodation and drowned_band at 0 Myr and at the start.
rem  Full-length cloud pre-check, seed 1022 (unused), control and 'accommodation', read only, not pooled (journal 129
rem   sec. 6; precheck-p96-spread\precheck_1022.txt): every veto reading passes: gates pooled young-ocean -118 / thin-crust
rem   -4 m (control -169 / +8); highest node 6,025 m (6,008); detached 0.0025 (0.0024); crust residual whole run +0.20
rem   km3/yr (+0.46), closure 0.0000. At -250 Myr partial r(S, b | Tc) -0.55 (+0.25); continents under the datum 0.133
rem   (0.258); continent at -1,000 ... -200 m 0.018 of the surface (0.046). At 0 Myr: land 0.291 (0.262; Earth 0.281);
rem   continent at -1,000 ... -200 m 0.047 (0.067); continents under the datum 0.234 (0.234), window 0.202 (0.216); the
rem   sea -118 m on the datum (-27); but land p50 575 m (528; Earth 380), window 614 (516). One seed; the sweep decides.
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
set A=start_cover-copula
set B=start_cover-accommodation
%PY% -u -m scripts.run_sweep --name start-cover --grid start_cover=copula,accommodation --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261120-190819\sediment_rock-burial || goto fail
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
%PY% %P96%\drowned_band.py --t=t-0250.0 "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\drowned_band_start.txt" 2>&1 || goto fail
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
