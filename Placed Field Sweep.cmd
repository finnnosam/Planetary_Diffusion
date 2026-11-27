@echo off
rem Placed Field Sweep (21 Nov; journal 132, ninety-ninth pass). Rule fixed 21 Nov before the runs.
rem  Why: STATE section 0, the continents' level on the datum; item 2, the dynamic field's placement (un-parked with the
rem  user, 21 Nov). Pass 98 (journal 131): without each body's own field the planet's continents (crystalline >= 30 km)
rem  stand at Earth's to ~40 m at every quantile (p5 ... p90 -325 / -142 / 134 / 413 / 764 / 1,340 m against -281 / -124 /
rem  138 / 454 / 775 / 1,391); the plate-blind measured field widens them on both sides (p10 -132, p90 +173 m; mean dyn on
rem  continents under / above the datum -448 / +111 m, Earth -136 / -13), and the drowned band's excess is normal crust
rem  inland on the field's lows (3.3% of the surface, dyn -443 m). Earth's field narrows Earth's continents.
rem  flow_power (built this pass; tests\test_flow_power.py): with dyn_rule 'slabs' (D, built 31 Oct), 'flow' keeps the
rem  slab flow's own amplitude (727 m rms, 3.3x Earth's whole observed signal: parked); 'measured' keeps the flow's
rem  pattern and rescales each degree to Earth's measured residual-topography power (Valentine & Davies 2020, l 1-30,
rem  air-loaded, exactly as dyn_rule 'measured' reads it). So the arm differs from the control in one thing: where the
rem  field's highs and lows lie (the planet's own slab history against random phases). No number added or fitted. The
rem  field is zero until the record first loads below 300 km (~7 Myr after the start). Default path identical step by
rem  step (n = 48, 70 steps, 71 of 71 states).
rem  Arms (2 x 6 seeds = 12 runs), dyn_rule, flow_power:
rem   A measured, flow    control (today's PARAMS; must equal sweep-20261121-095209's marine_diffusion-none)
rem   B slabs, measured   the measured field placed by the slab record
rem  Rule. Vetoes on B: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close: on every seed
rem   the series' crust ledger closes (v3_split closure <= 0.001 km3/yr) and crust_budget's crust residual is the
rem   control's size by seed; a seed outside the control's range passes if v3_split puts the excess in merges / pieces /
rem   holes, or its update part, as a share of crystalline erosion, lies within this sweep's control seeds' range widened
rem   by 2 sd of them on each side; V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the
rem   highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: B is adopted if it passes every veto (a realism change: the plate-blind random phases are replaced by the
rem   planet's own slab history at the same measured power per degree; no number added). Fails a veto: the control
rem   kept, the failing veto is next.
rem  Design check (read, not a veto): dynamic_topography's P-b, corr(field, distance to young slabs) > 0 on the
rem   six-seed mean (lows over slabs), and the air-loaded rms at Earth's 0.41 km (by construction).
rem  Expected (not a band): the continents' spread narrows (dyn_tail: z against z - dyn); the prediction carried from
rem   journal 74: placement may lower the continents' level (Earth's continents sit on its lows, -69 m on Holdt's field),
rem   so read the sea on the datum, the young floor (basin) and land beside the spread; a lower level is scored, the next
rem   item, not a reason to reject. Caveat: the arm puts the slab pattern at every degree; Earth's short wavelengths
rem   (l >~ 15) are mostly lithosphere and upper mantle (Davies et al. 2019), which the planet carries as C.
rem  Scored, reported (no veto): as the Marine Diffusion Sweep (level_split, rest_split, continent, three_factor, basin,
rem   hyps, column_split per arm; coast_profile / coast_budget, lith_crust2, thin_where, regime_split, sediment_where,
rem   sediment_by_height, thick_tail, crust_budget, v3_per_seed, v3_events, v3_split, the margins band, passive_survival,
rem   cover_accommodation, drowned_band, shelf_profile, drowned_setting), plus dyn_tail, dyn_inland, drowned_split,
rem   drowned_land and dynamic_topography.
rem  Full-length cloud pre-check, seed 1024 (unused), control and B, read only, not pooled (journal 132 sec. 4;
rem   precheck-p99-field\precheck_1024.txt): every veto reading passes on B
rem   (the second build; the first failed V5 at the start, journal 132 sec. 3): gates pooled young-ocean -76 / thin-crust
rem   +24 m (control -120 / +13); highest node in the snapshots 5,949 m (6,273); detached 0.0025 (0.0018); crust residual
rem   whole run / last 40 Myr -0.04 / +0.11 km3/yr (+0.06 / +0.48); the start identical to the control's. At 0 Myr:
rem   continents (crust >= 30 km) under the datum 0.180 (0.317), p10-p90 on the datum -78 / 89 / 418 / 1,114 / 1,960 m
rem   (-568 / -65 / 206 / 622 / 1,270); drowned band 0.036 (0.058); land 0.273 (0.268), land p50 / p90 509 / 2,009 m
rem   (385 / 1,417); sea on the datum -55 m (-92); young floor 0.382 (0.367); field rms 416 m, corr(field, distance to young
rem   slabs) +0.46; plates at 0 Myr 47 (37). The low tail rises, the high side rises too on this seed. One seed; the
rem   sweep decides.
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
set P98=%CAL%\precheck-p98-low
set A=dyn_rule-measured_flow_power-flow
set B=dyn_rule-slabs_flow_power-measured
%PY% -u -m scripts.run_sweep --name placed-field --grid dyn_rule,flow_power=measured:flow,slabs:measured --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261121-095209\marine_diffusion-none || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% -u -m planet_reviews.dynamic_topography "%SW%\%A%" "%SW%\%B%" --json "%SW%\reviews\dynamic_topography.json" > "%SW%\reviews\dynamic_topography.txt" 2>&1 || goto fail
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
%PY% %P96%\dyn_tail.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\dyn_tail.txt" 2>&1 || goto fail
%PY% %P96%\dyn_inland.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\dyn_inland.txt" 2>&1 || goto fail
%PY% %P97%\shelf_profile.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\shelf_profile.txt" 2>&1 || goto fail
%PY% %P97%\drowned_setting.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\drowned_setting.txt" 2>&1 || goto fail
%PY% %P98%\drowned_split.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\drowned_split.txt" 2>&1 || goto fail
%PY% %P98%\drowned_land.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\drowned_land.txt" 2>&1 || goto fail
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
