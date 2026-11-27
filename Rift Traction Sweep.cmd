@echo off
rem Rift Traction Sweep (21 Nov; journal 139, pass 103). Rule fixed 21 Nov before the runs.
rem  Why: STATE section 0, the continents' level on the datum; item 3, how the dynamic field meets the continents
rem  (journal 137 sec. 5). Pass 103's probes (outputs\planet-calibration\precheck-p103-meet):
rem  - earth_cont_corr.py: Earth's residual topography (VD2020 / Holdt 2022) against Earth's continents (CRUST1.0), band
rem    correlation l 1-3 -0.62 / -0.31 (over oceanic nodes alone -0.69 / -0.42: the oceans constrain it; l 2 -0.76 /
rem    -0.83), l 4-10 -0.06 / -0.09 (oceans -0.06 / +0.01), l 11-30 ~0. Earth's field knows the continents at l 1-3 only.
rem  - planet_cont_corr.py (six seeds x 0 / -20 / -40 Myr): today's model l 1-3 / 4-10 / 11-30 +0.27 / -0.04 / -0.03
rem    (l 1 +0.69): the l >= 4 part already matches Earth's zero, so imposing Earth's correlation on the plate-blind rest
rem    (candidate (a)) has nothing to impose; at l 1-3, the slab record's degrees, the sign is the opposite of Earth's.
rem  - cont_spectrum.py: every model since pass 98 keeps its continents in one cap: degree-1 share of the continent
rem    function 0.52-0.54 and resultant 0.43-0.44 at 0 Myr (start 0.71 / 0.59; Earth 0.25 / 0.31); continent on one
rem    plate 65-104 Mkm2 (Earth 42). The cap sits on its slab girdle's interior high.
rem  - traction_read.py: under the plate holding the most continent, the placed flow (scaled like its topography to
rem    Earth's power) diverges on 5 of 6 seeds, +0.35 ... +1.17 cm/yr outward; its traction through the existing basal
rem    and keel drag is 0.5-4 x 10^20 N, one to eight 1,000 km trenches' slab pull. The split rule's forces leave it out.
rem  - cloud split logs (seed 1001, both arms of the last sweep): a plate-wide cut is accepted with p ~0.01-0.05 per
rem    step (1 / L), a corner with 0.03-0.36; an unbroken supercontinent is a draw of that lottery (~8%).
rem  rift_traction (built this pass; tests\test_rift_traction.py): False (today) / True: the split rule's forces
rem  (rift_line 'force': the ranking, the tension test, the cut chosen) add the basal traction of the planet's own mantle
rem  flow, C_b (1 + (keel - 1) c) v_flow, the flow scaled with its topography (MantleFlow.placed_velocity). Split rate and
rem  plate motions unchanged. No number added. Default path identical step by step at threads 1 (71 of 71 states);
rem  threads 2 cannot be checked in the cloud (HEAD against HEAD differs there); the borrow's check seed tests it here.
rem  Arms (borrowed control: 1 + 6 runs), rift_traction:
rem   A False  control (today's PARAMS); its first seed re-run and compared with sweep-20261121-185217's
rem            flow_placed-earth_corr, the other five hard-linked from it
rem   B True   the planet's own flow pulls on the plate in the split rule's forces
rem  Rule. Vetoes on B: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close: on every seed
rem   the series' crust ledger closes (v3_split closure <= 0.001 km3/yr) and crust_budget's crust residual is the
rem   control's size by seed; a seed outside the control's range passes if v3_split puts the excess in merges / pieces /
rem   holes, or its update part, as a share of crystalline erosion, lies within this sweep's control seeds' range widened
rem   by 2 sd of them on each side; V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the
rem   highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: B is adopted if it passes every veto (a realism change: a force the model computes and Earth's
rem   supercontinents feel enters the forces that crack its plates; no number added). Fails a veto: the control kept,
rem   the failing veto is next.
rem  Design check (read, not a veto): stats.json chord_stats (the share of 'force' cuts) higher on B; traction_read's
rem   outward flow under the biggest continental plate.
rem  Expected (not a band): the continent on one plate (plate_continent), the degree-1 share and the resultant
rem   (cont_spectrum) fall toward Earth's 42 Mkm2 / 0.25 / 0.31; planet_cont_corr's l 1-3 toward Earth's negative; the
rem   interiors' field (degree_split, dyn_tail) down. Risk, scored: more and larger young oceans (basin's young floor
rem   0.391 against 0.30), T2-T5, the rifted production (per unit continent, 1.13x; the rule keeps it by design).
rem  Scored, reported (no veto): as the Placed Share Sweep, plus drowned_far, plate_continent, planet_cont_corr,
rem   cont_spectrum and traction_read.
rem  Full-length cloud pre-check, seed 1027 (unused), control and B, read only, not pooled (journal 139 sec. 5;
rem   precheck-p103-meet\precheck_1027.txt, precheck_1027_breakup.txt, crust_budget_1027.txt): every veto reading passes
rem   on B (gates -78 / +23 m, highest node 6,181 m, detached 0.0030, crust residual -0.03 / +0.01 km3/yr). From -250 to
rem   -150 Myr every draw on the supercontinent plate finds a cut in tension (control 0.66-0.78), its median part 75-78
rem   Mkm2 (22-38); it splits at -200 ... -150 Myr. At 0 Myr: continent on the three largest plates 49 / 42 / 39 Mkm2
rem   (control 117 / 28 / 21; Earth 42 / 35 / 27); degree-1 share 0.30 (0.53; Earth 0.25); resultant 0.31 (0.44; 0.31).
rem   Costs on this seed: land 0.259 (0.306), continents under the datum 0.210 (0.193), p10 -389 m (-154), submerged
rem   continent 0.148 (0.101), young floor 0.400 (0.342); planet_cont_corr l 1-3 +0.23 (+0.20). One seed; the sweep decides.
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
set P102=%CAL%\precheck-p102-share
set P103=%CAL%\precheck-p103-meet
set A=rift_traction-False
set B=rift_traction-True
%PY% -u -m scripts.run_sweep --name rift-traction --grid rift_traction=False,True --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261121-185217\flow_placed-earth_corr --borrow-control || goto fail
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
%PY% %P102%\drowned_far.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\drowned_far.txt" 2>&1 || goto fail
%PY% %P102%\plate_continent.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\plate_continent.txt" 2>&1 || goto fail
%PY% %P103%\planet_cont_corr.py "%SW%\%A%" "%SW%\%B%" --tags t+0000.0,t-0020.0,t-0040.0 > "%SW%\reviews\planet_cont_corr.txt" 2>&1 || goto fail
%PY% %P103%\cont_spectrum.py "%SW%\%A%" "%SW%\%B%" > "%SW%\reviews\cont_spectrum.txt" 2>&1 || goto fail
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
%PY% %P103%\traction_read.py "%SW%\%ARM%" > "%SW%\reviews\traction_read_%1.txt" 2>&1 || exit /b 1
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
