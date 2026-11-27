@echo off
rem Rock Cover Sweep (20 Nov; journal 126, ninety-fifth pass). Rule fixed 20 Nov before the runs.
rem  Why: STATE section 0, the continents' level on the datum. On today's model (water_load_rule 'column', adopted this
rem  pass from sweep-20261119-215825) pass 94's rule (precheck-p94-level\level_split.py, fixed before its runs) ranks one
rem  candidate over 30-45 km crust: (b) the sediment column -195 m (2 sem 47); C -65 and dyn -6 do not rank. Its sign is
rem  the cover's place: the planet carries 4.6 / 3.5 km of sediment on 20-25 / 25-30 km crust (Earth 2.6 / 2.8) and
rem  0.75 / 0.49 / 0.32 km on 35-40 / 40-45 / 45-55 km (Earth 1.32 / 1.62 / 1.08). Mechanism: under
rem  sediment_entrainment 'cover' the whole sediment layer is loose (entrained at theta 0 down to its base), also the
rem  start's continental cover, which is Earth's stable platform cover: Palaeozoic-Mesozoic sedimentary rock. On the
rem  planet it leaves the high land within ~40 Myr and settles low (journal 116). On Earth sediment buried below ~2 km
rem  is rock (quartz cementation from 70-80 C, Walderhaug 1994; mechanical compaction to 2-3 km, Bjorlykke 1999).
rem  sediment_rock (built this pass; tests\test_sediment_rock.py): 'none' as built; 'burial': each plate carries a rock
rem  layer (solid m, <= its sediment); sediment buried below rock_depth_m 2,000 m (bulk) turns to rock; the start's
rem  continental cover is rock; only the loose top is the 'cover' law's cover, rock erodes by the bedrock law (theta
rem  and K as crystalline: no new erosion number). Default path identical step by step (n = 48, 70 steps, threads 1).
rem  Arms (2 x 6 seeds = 12 runs), sediment_rock:
rem   none    control (today's PARAMS; must equal sweep-20261119-215825's water_load_rule-column)
rem   burial  sediment buried below 2 km and the start's cover are rock
rem  Rule. Vetoes on 'burial': V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close: on
rem   every seed the series' crust ledger closes (v3_split closure <= 0.001 km3/yr; the series carries sim.ledger since
rem   this pass) and crust_budget's crust residual is the control's size by seed; a seed outside the control's range
rem   passes if v3_split puts the excess in merges / pieces / holes, or its update part, as a share of crystalline
rem   erosion, lies within this sweep's control seeds' range widened by 2 sd of them on each side (the update's net
rem   varies in sign; the Water Load Sweep's column arm: -0.016 ... +0.080 outside merge steps, sd 0.035); V4
rem   plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the highest node over the run
rem   <= 8,850 m or no higher than the control's.
rem  Adoption: 'burial' is adopted if it passes every veto (a mechanism with a measured onset and no fitted value:
rem   Earth's old cover is rock, its fresh sediment loose). Fails a veto: the control kept, the failing veto is next.
rem  Design check (read, not a veto): level_split's weighted sediment gap over 30-45 km nearer zero than the control's
rem   (-194 m on the Water Load Sweep's column arm), and the mean sediment on 35-40 and 40-45 km crust higher.
rem  Expected (not a band): continents on the datum up; land up; the sea down; risk (the pre-check): the high land
rem   keeps its cover and stands above Earth's (land p50, the spread), crystalline erosion up, crust alone past -3%.
rem  Scored, reported (no veto): level_split and rest_split per arm (precheck-p94-level), continent, three_factor,
rem   basin, hyps, column_split per arm, coast_profile / coast_budget, lith_crust2, thin_where, regime_split,
rem   sediment_where, sediment_by_height, thick_tail, crust_budget, v3_per_seed, v3_events, v3_split (precheck-p95-v3), the
rem   margins band, passive_survival.
rem  Full-length cloud pre-check, seed 1020 (unused), control and 'burial', read only, not pooled (journal 126 sec. 6):
rem   every veto reading passes on 'burial': gates pooled young-ocean -117 / thin-crust +35 m (control -99 / +34);
rem   highest node 6,243 m (control 6,243, both at the start); detached <= 0.0062 (0.0043); crust residual whole run
rem   +0.26 km3/yr (control +0.04): outside merge steps 8.7% of crystalline erosion (control 6.5%). Continents on the
rem   datum p25 / p50 / p75 79 / 470 / 1,056 m, under the datum 0.187 (control 59 / 258 / 658, 0.218); land 0.262
rem   (0.236); the sea +51 m on the datum (+109); drowned flagged 0.34 (0.44). But land p50 604 m (331; Earth 380):
rem   the high land keeps its cover; crystalline erosion -3.27 km3/yr (-2.91), crust alone -3.2% (-2.7%). The cover by
rem   bin moves little toward Earth's (35-40 / 40-45 / 45-55 km 827 / 559 / 462 m, control 653 / 462 / 268, Earth
rem   1,317 / 1,616 / 1,076) and the thin bins hold more (5.5 / 4.9 km on 20-25 / 25-30). One seed; the sweep decides.
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
set A=sediment_rock-none
set B=sediment_rock-burial
%PY% -u -m scripts.run_sweep --name rock-cover --grid sediment_rock=none,burial --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261119-215825\water_load_rule-column || goto fail
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
