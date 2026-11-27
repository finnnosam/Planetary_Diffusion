@echo off
rem Placed Share Sweep (21 Nov; journal 136, pass 101). Rule fixed 21 Nov before the runs.
rem  Why: STATE section 0, the continents' level on the datum; item 3, the placed field's degree content (journal 134
rem  sec. 5). Pass 101's probes (outputs\planet-calibration\precheck-p101-degree):
rem  - degree_split.py on today's model (sweep-20261121-133838, arm B): the drowned normal crust's low is the field's
rem    l 4-10 (-249 m of -293), the interiors' high is l 1-3 and 4-10 (+114 / +76 m at 600-1200 km, +143 / +145 beyond);
rem    l 11-31 is near zero on every set. Earth (Holdt 2022): interiors +139 / -167 m at l 1-10 / 11-30 (beyond 1,200 km
rem    -8 / -354): Earth's low-degree field is not low over the continents' margins and its interior lows are short.
rem  - trench_field.py: continent (crystalline >= 30 km) within 0-300 / 300-600 / 600-1000 km of subduction sits on
rem    l 1-10 field -376 / -156 / -91 m on the planet against Earth's +36 / +55 / +34 m.
rem  - earth_degree_corr.py: the degree correlation of Earth's measured residual topography (Valentine & Davies 2020)
rem    with the same flow on Earth's own subduction history (E3: Muller et al. 2019 through D's record and kernels, journal
rem    70, slab_alpha thermo): l 1 / 2 / 3 r +1.00 / +0.88 / +0.32; bands l 4-10 +0.07, 11-30 +0.00 (the tomographic
rem    slab proxy E1: +0.21 / +0.09). On Earth the slabs account for degrees 1-2 and little else; the planet's flow_power
rem    'measured' gives the slab pattern all of Earth's power at every degree.
rem  - Caveat (earth_pb.py, planet_pb.py): corr(field, distance to subduction) by band, Earth (VD2020 mean, Bird's SUB)
rem    l 1-3 / 4-10 / 11-30 +0.46 / +0.23 / -0.00 (Holdt 1-10 / 1-30 +0.15 / +0.17), the planet today +0.32 / +0.21 / +0.05:
rem    Earth's l 4-10 does have some lows near trenches, which E3's slab flow does not reproduce (r +0.07); under B the
rem    planet's l 4-10 loses its trench correlation. Read on both arms; the arm rests on the degree correlation.
rem  flow_placed (built this pass; tests\test_flow_placed.py): 'all' (today) places the whole measured power per degree;
rem  'earth_corr' places max(r_l, 0)^2 of it (Earth's measured share: 0.99 / 0.78 / 0.10 at l 1-3, < 0.05 of the power
rem  at l 4-30) and draws the rest, (1 - share), as a plate-blind field like dyn_rule 'measured''s (same stream, after the
rem  start's field; blended between knots at constant variance). Total power per degree stays Earth's. The start is the
rem  control's. No number fitted. Default path identical step by step (n = 48, 70 steps, 71 of 71 states).
rem  Arms (borrowed control: 1 + 6 runs), flow_placed:
rem   A all         control (today's PARAMS); its first seed re-run and compared with sweep-20261121-133838's
rem                 dyn_rule-slabs_flow_power-measured, the other five hard-linked from it
rem   B earth_corr  the slab record places Earth's measured share of each degree; the rest plate-blind
rem  Rule. Vetoes on B: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close: on every seed
rem   the series' crust ledger closes (v3_split closure <= 0.001 km3/yr) and crust_budget's crust residual is the
rem   control's size by seed; a seed outside the control's range passes if v3_split puts the excess in merges / pieces /
rem   holes, or its update part, as a share of crystalline erosion, lies within this sweep's control seeds' range widened
rem   by 2 sd of them on each side; V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the
rem   highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: B is adopted if it passes every veto (a realism change: the share of the field the slabs place becomes
rem   the share Earth's own slabs account for, measured; no number added). Fails a veto: the control kept, the failing
rem   veto is next.
rem  Design check (read, not a veto): dynamic_topography's air-loaded rms at Earth's 0.41 km (by construction); P-b on
rem   the whole field > 0 and below A's +0.42; P-b of the placed part alone > 0; planet_pb by band beside Earth's.
rem  Expected (not a band): degree_split's l 4-10 low on the drowned normal crust and its highs on the interiors go
rem   toward zero on the six-seed mean; trench_field's l 1-10 near subduction toward Earth's. Risk, scored: the plate-blind
rem   rest widens the continents as the old plate-blind default did (p10 -326 m against B's -214), so read continent,
rem   three_factor, dyn_tail and hyps beside the drowned band; a lower p10 is scored, the next item, not a reason to reject.
rem  Scored, reported (no veto): as the Placed Field Sweep (level_split, rest_split, continent, three_factor, basin, hyps,
rem   column_split per arm; coast_profile / coast_budget, lith_crust2, thin_where, regime_split, sediment_where,
rem   sediment_by_height, thick_tail, crust_budget, v3_per_seed, v3_events, v3_split, the margins band, passive_survival,
rem   cover_accommodation, drowned_band, shelf_profile, drowned_setting, dyn_tail, dyn_inland, drowned_split, drowned_land,
rem   dynamic_topography), plus degree_split, trench_field and planet_pb.
rem  Full-length cloud pre-check, seed 1026 (unused), control and B, read only, not pooled (journal 136 sec. 4;
rem   precheck-p101-degree\precheck_1026.txt): every veto reading passes on B (gates -94 / +6 m, highest node the start's
rem   6,457 m, detached 0.0000, crust residual +0.16 / +0.26 km3/yr). At 0 Myr B is worse on this seed: continents under
rem   the datum 0.324 (control 0.227), p10 -699 m (-297); its plate-blind l 4-10 puts -225 m under the drowned crust. One
rem   seed; the sweep decides.
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
set P101=%CAL%\precheck-p101-degree
set A=flow_placed-all
set B=flow_placed-earth_corr
%PY% -u -m scripts.run_sweep --name placed-share --grid flow_placed=all,earth_corr --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261121-133838\dyn_rule-slabs_flow_power-measured --borrow-control || goto fail
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
%PY% %P101%\degree_split.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\degree_split.txt" 2>&1 || goto fail
%PY% %P101%\trench_field.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\trench_field.txt" 2>&1 || goto fail
%PY% %P101%\planet_pb.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\planet_pb.txt" 2>&1 || goto fail
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
