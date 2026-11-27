@echo off
rem Collision Lithosphere Sweep (8 Nov; journal 101, seventy-seventh pass). Rule fixed 8 Nov before the runs.
rem  Why: the land's height is the largest failure (land p50 ~1,200 m against 384) and half of it is the column
rem  (journal 98). Lamb et al. 2020 (data\Lamb2020\ds02.xyz: continental means of measured lithosphere, Priestley and
rem  McKenzie 2016, on CRUST1.0 crust) show why Earth's land is flat against crust thickness: its lithosphere thickens
rem  with its crust, 6.7-7.4 km per km (thin-crust regions 80-92 km, thick-crust interiors 170-192 km, Eurasia with
rem  Tibet 188), and once reduced for it Earth's elevation rises 128 m per km of crust, as the planet's crust column.
rem  The planet's lithosphere does not follow its crust (corr 0.04, 0.24 km per km; Deep Crust Sweep, adopted arm,
rem  lith_crust.txt): C resets it to age 0 (L 93.6 km) at collision contacts, so crust of 45-80 km sits on L p50 105 km.
rem  collision_lithosphere 'thicken': pure shear, the mirror of C's stretching rule: where collision thickens the crust
rem  by T_after / T_before the lithosphere thickens by the same factor, with no reset (the continental arc belt still
rem  resets). A mechanism; Lamb's 30 kg/m3 and Artemieva's L(t) unchanged. No new number.
rem  Second run (8 Nov, journal 102). The first (sweep-20261108-190402) stopped: on the thicken arm L ratcheted without
rem  bound (collision volume each step while the collapse flow and foundering held the crust near 50-65 km; L > 300 km on
rem  1% of the continent by -40 Myr, then ~1e10 km, elevations and the sea level non-finite on 4 seeds). Fixed by two
rem  mechanisms Earth has: lithosphere thickened past the thickest stable root is removed (convective removal, Houseman
rem  et al. 1981; lith_max_km 237.6 = C's own law at 3.6 Ga), and a foundering root (eclogite_rule) takes its mantle
rem  lithosphere (delamination: age reset). Cloud pre-check, seed 1007, 250 Myr (read only, not pooled): L max over the
rem  run 243 km (bound + aging); land p50 1,095 -> 1,021 m, p90 2,534 -> 2,191; dL/dT 0.18 -> 0.28 km per km (thickening
rem  at collisions does not make Earth's L-T relation: the thickest crust delaminates, normal crust is untouched).
rem  Rule unchanged.
rem  Expected (written before the runs, not a band): thick crust on thicker lithosphere (L p50 on 45-80 km crust from
rem  ~105 toward 150-200 km), its land 0.3-0.8 km lower, land p50 down ~100-250 m, dL/dT up from 0.2 toward 2-4 km
rem  per km (pure shear alone gives L / T; Earth's 7 includes cratons formed thick: the start's join is the next
rem  candidate if this falls short), the column share (column_split) down by a similar amount. Thicker belts are
rem  stronger under rift_line's strength map: rifting may move.
rem  Arms (2 x 6 seeds = 12 runs, about 40 min); folders collision_lithosphere-...:
rem   reset     today's PARAMS (control; must equal crust_density_rule-pressure_crust_kappa_m2_yr-1000 of
rem             sweep-20261108-170150: identity.txt)
rem   thicken   the candidate
rem  Rule. Vetoes, every one on its arm: V1 pytest passes (incl. test_collision_lithosphere); V2' elevation gates
rem   +-300 m on 6 of 6 seeds; V3 ledgers close; V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on
rem   every seed the highest node over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: collision_lithosphere 'thicken' if its arm passes every veto (pure shear, as the stretching rule
rem   already does; it moves the planet toward Earth's measured lithosphere-crust relation).
rem  Scored, reported (no veto): lith_crust (dL/dT, corr, L by crust bin, Lamb-reduced elevation) and column_split per
rem   arm, land fraction and p25 / p50 / p90, sea, C1, crust bands, export, volume, K1-K4, edge, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P77=%CAL%\precheck-p77-lithosphere
%PY% -u -m scripts.run_sweep --name collision-lithosphere --grid collision_lithosphere=reset,thicken --seeds 1001,1002,1003,1004,1005,1006 --control collision_lithosphere-reset --identity %CAL%\sweep-20261108-170150\crust_density_rule-pressure_crust_kappa_m2_yr-1000 || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P77%\lith_crust.py "%SW%\collision_lithosphere-reset" "%SW%\collision_lithosphere-thicken" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P75%\column_split.py "%SW%\collision_lithosphere-reset" > "%SW%\reviews\column_split_reset.txt" 2>&1 || goto fail
%PY% %P75%\column_split.py "%SW%\collision_lithosphere-thicken" > "%SW%\reviews\column_split_thicken.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
