@echo off
rem Deep Crust Sweep (8 Nov; journal 100, seventy-sixth pass). Rule fixed 8 Nov before the runs.
rem  Why (journal 98, the step-back; order chosen by Claude, the user's delegation of 8 Nov): the land's height is the
rem  largest failure (p50 1,170 m against 384) and half of it is the column: at a given crust thickness the planet's
rem  land stands too high, rising ~120 m per km of crust from 30 to 45 km against Earth's ~25. This sweep tests the
rem  first mechanism for it: the model's crust has one density at every depth (rho_crust 2,830), so a 60 km root is
rem  as buoyant per km as 35 km crust right up to where eclogite_rule founders it (65 km, where mafic crust reaches
rem  the mantle's density, Chapman et al. 2019). crust_density_rule 'pressure' (planet_tectonics.deep_crust): the
rem  deepest 0.324 of the crystalline crust (CRUST1.0 types, 2,947 kg/m3; measured) gains density linearly in pressure
rem  from 1.0 GPa (garnet in, Green and Ringwood 1967) to the mantle's 3,300 at 2.0 GPa (Chapman et al. 2019), held
rem  beyond; a dry flexural load. Equilibrium assumed (the upper bound). No fitted number.
rem  Second factor, read only: the collapse flow (crust_kappa_m2_yr 1,000, Tibet-tuned, fitted) off, as in journal 77
rem  (then: crust > 45 km 0.106, rejected on export and volume when those were vetoes; not swept since the 4 Nov rules).
rem  Expected (written before the runs, not a band; static estimate on the pelagic arm's 0 Myr surfaces): 'pressure'
rem  lowers land p50 ~30 m and p90 ~240 m, the 45-55 km bin ~250 m and the 55-80 km bin ~800 m (away from Earth's
rem  4,356 m: scored), the column share (column_split) by tens of metres only: the interior's share needs another
rem  mechanism (journal 100). Collapse off: less thick crust inland, steeper belts, more export; towers possible (V5).
rem  Arms (4 x 6 seeds = 24 runs, about 80 min); folders crust_density_rule-..._crust_kappa_m2_yr-...:
rem   uniform / 1000    today's PARAMS (control; must equal pelagic_m_myr-1.1 of sweep-20261108-150308: identity.txt)
rem   pressure / 1000   the candidate
rem   uniform / 0       collapse flow off (read only)
rem   pressure / 0      both (read only)
rem  Rule. Vetoes, every one on its arm: V1 pytest passes (incl. test_deep_crust, test_pelagic); V2' elevation gates
rem   +-300 m on 6 of 6 seeds (elevation_check reads the new load); V3 ledgers close; V4 plate_pieces_review T1 detached
rem   area <= 0.01 on every seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the
rem   control's there.
rem  Adoption: crust_density_rule 'pressure' if the pressure / 1000 arm passes every veto (a process Earth has, at
rem   measured values; it makes the root's density consistent with the foundering rule). The collapse-off arms are not
rem   adopted here whatever they show; their reading orders the thick-tail item.
rem  Scored, reported (no veto): column_split per arm (land with the sea at the datum on Earth's thickness weights; z by
rem   thickness bin on the datum), land fraction and p25 / p50 / p90, sea, crust bands, export, volume, coast_budget (C1),
rem   K1-K4, edge, margins band, T2-T5, plates, highest node.
rem  Runs with scripts.run_sweep: tests, the 24 runs and the reviews under one progress line; then column_split per arm
rem  into <sweep>\reviews\.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
%PY% -u -m scripts.run_sweep --name deep-crust --grid crust_density_rule,crust_kappa_m2_yr=uniform:1000,pressure:1000,uniform:0,pressure:0 --seeds 1001,1002,1003,1004,1005,1006 --control crust_density_rule-uniform_crust_kappa_m2_yr-1000 --identity %CAL%\sweep-20261108-150308\pelagic_m_myr-1.1 || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
for %%a in (uniform_crust_kappa_m2_yr-1000 pressure_crust_kappa_m2_yr-1000 uniform_crust_kappa_m2_yr-0 pressure_crust_kappa_m2_yr-0) do (
  %PY% %P75%\column_split.py "%SW%\crust_density_rule-%%a" > "%SW%\reviews\column_split_%%a.txt" 2>&1 || goto fail
)
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
