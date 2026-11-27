@echo off
rem Thick Tail Sweep (14 Nov; journal 109, eighty-fourth pass). Rule fixed 14 Nov before the runs.
rem  Why: STATE section 0, item 3, the thick tail. The planet keeps 0.174 of its continent above 45 km (Earth 0.086); its
rem  plateaus erode slowly, keep their sediment high (the sea's supply), stand high (the continents' spread) and hold
rem  continental volume off the area (journal 108). The crust flow that shapes them is fitted: kappa_normal 30 plus a
rem  collapse term above 45 km whose four numbers were tuned to Tibet (crust_kappa_m2_yr 1000, collapse_thick_m 45 km,
rem  collapse_scale_m 6 km, collapse_max_factor 200). crust_flow_rule replaces those four by the flow the crust's own
rem  rheology gives (planet_tectonics.crust_flow): a steady geotherm from the column's crust T and C's lithosphere L
rem  (heat production 1.0 / 0.4 uW/m3 in the felsic / lower crust, Hasterok & Chapman 2011; the LAB at 1,175 C, the
rem  isotherm Hoggard et al. 2020 map as L), dislocation creep of wet quartzite (Gleason & Tullis 1995) and wet anorthite
rem  in the lower 0.324 (Rybacki & Dresen 2000), driven by the isostatic surface slope; no number fitted to the planet.
rem  Arms (4 x 6 seeds = 24 runs); folders crust_flow_rule-...:
rem   fitted   today's PARAMS (control; must equal margin_wedge-True of sweep-20261114-104046: identity.txt)
rem   channel  fixed boundary everywhere (zero velocity at the surface and the Moho; Rey et al. 2001's convergent mode)
rem   sheet    free boundary everywhere (the crust spreads from its surface down; the post-convergent mode)
rem   rey      each mode where Rey et al. put it: fixed within 290 km of active convergence (arc belt, closing collision
rem            boundary; half of Earth's collisional thickening lies within 290 km, crust_kernels.COLLISION_KERNEL), free
rem            elsewhere
rem  Rule. Vetoes on each candidate arm: V1 pytest passes (incl. test_crust_flow); V2' elevation gates +-300 m on 6 of 6
rem   seeds; V3 ledgers close (crust_budget residuals of the control's size); V4 plate_pieces_review T1 detached area <=
rem   0.01 on every seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption, in this order (by mechanism, not by score): rey if it passes every veto (both modes where they apply);
rem   else sheet (the planet's thick crust is post-collisional: 98% collided, kept for 100s of Myr, journal 75); else
rem   channel. None passing: fitted kept, and the failing veto is the next thing to work on.
rem  Full-length cloud pre-check, seed 1015 (unused), fitted / channel / sheet / rey, read only, not pooled
rem   (precheck-p84-thicktail\precheck_read_1015.txt, thick_tail_1015.txt; journal 109 section 5): bounded on all four,
rem   no towers (highest 6,589 / 6,471 / 5,490 / 6,944 m), gates young ocean -96 / -137 / -75 / -231 m, thin crust
rem   -110 / -106 / -97 / -172 m; > 45 km 0.199 / 0.185 / 0.258 / 0.259; volume -1.2 / -1.3 / -2.7 / -3.8%; to the floor
rem   0.66 / 0.68 / 0.70 / 1.41 km3/yr. The power law flattens plateaus at 45-55 km instead of thinning them.
rem  Expected (written before the sweep's runs, not a band): sheet and rey flatten the tail without thinning it (> 45 km
rem   up, > 65 km gone under sheet); rey sends more sediment to the floor at a volume cost, its young-ocean gate the
rem   nearest the band; channel near the control. Land p50 unresolved. The new arms take ~1.5x a control run's time.
rem  Scored, reported (no veto): thick_tail (crust bands, > 45 / > 55, Moho temperature, the flow's substeps and
rem   ceiling share), interior_probe (slope of land > 2 km, depth of > 45 km crust from its edge), three_factor, basin,
rem   continent, hyps, column_split, sediment_where, sediment_by_height per arm; lith_crust2; land fraction and
rem   p25 / p50 / p90, C1, crust bands, export, volume, K1-K4, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
set P83=%CAL%\precheck-p83-sediment
set P84=%CAL%\precheck-p84-thicktail
%PY% -u -m scripts.run_sweep --name thick-tail --grid crust_flow_rule=fitted,channel,sheet,rey --seeds 1001,1002,1003,1004,1005,1006 --control crust_flow_rule-fitted --identity %CAL%\sweep-20261114-104046\margin_wedge-True || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P80%\lith_crust2.py "%SW%\crust_flow_rule-fitted" "%SW%\crust_flow_rule-channel" "%SW%\crust_flow_rule-sheet" "%SW%\crust_flow_rule-rey" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\crust_flow_rule-fitted" "%SW%\crust_flow_rule-channel" "%SW%\crust_flow_rule-sheet" "%SW%\crust_flow_rule-rey" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
for %%a in (fitted channel sheet rey) do (
  %PY% %P82%\three_factor.py "%SW%\crust_flow_rule-%%a" > "%SW%\reviews\three_factor_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\basin.py "%SW%\crust_flow_rule-%%a" > "%SW%\reviews\basin_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\continent.py "%SW%\crust_flow_rule-%%a" > "%SW%\reviews\continent_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\hyps.py "%SW%\crust_flow_rule-%%a" > "%SW%\reviews\hyps_%%a.txt" 2>&1 || goto fail
  %PY% %P75%\column_split.py "%SW%\crust_flow_rule-%%a" > "%SW%\reviews\column_split_%%a.txt" 2>&1 || goto fail
)
%PY% %P83%\sediment_where.py "%SW%\crust_flow_rule-fitted" "%SW%\crust_flow_rule-channel" "%SW%\crust_flow_rule-sheet" "%SW%\crust_flow_rule-rey" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\crust_flow_rule-fitted" "%SW%\crust_flow_rule-channel" "%SW%\crust_flow_rule-sheet" "%SW%\crust_flow_rule-rey" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
