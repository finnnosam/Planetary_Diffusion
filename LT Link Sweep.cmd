@echo off
rem L-T Link Sweep (11 Nov; journal 105, eighty-first pass). Rule fixed 11 Nov before the runs.
rem  Why: the land's height is the largest failure (STATE section 0); the column is its largest share and its cause is
rem  the lithosphere following the crust too little (Earth's SL2013sv L by crust bin on the nodes 97 / 98 / 107 / 143 /
rem  158 / 157 / 197 km, 20-25 ... 55-80 km; the lamb arm of sweep-20261109-222925 at 0 Myr 107 / 122 / 140 / 147 /
rem  145 / 139 / 155). Measured this pass (precheck-p81-ltlink): 82% of the 0 Myr continent is still on start
rem  lithosphere, whose L is flat across crust (146-154 km on 25-55 km); the young lithosphere under crust > 45 km is
rem  97 km (24% of the thick crust; half of it reset by the arc within the last Myr while margin shortening thickens it).
rem  Earth: L follows crust mostly through distance from the edge; within a distance bin rho 0.05-0.37 (0.10 beyond
rem  1,000 km, where 75% of the supercontinent lies); thick young crust on thick L (Andes 166 / 202, Tibet 158 / 196,
rem  the Cordillera's back-arc 119 on 45-55 km crust), thin arcs on thin L (72-99 km).
rem  Arms (2 x 2 x 6 seeds = 24 runs, about 80 min); folders start_lithosphere-..._belt_lithosphere-...:
rem   distance / reset   today's PARAMS (control; must equal lith_frame-lamb of sweep-20261109-222925: identity.txt)
rem   crust / reset      the start's ages coupled to its own crust within each distance bin (Gaussian copula, Earth's
rem                      rho per bin, lithosphere.EARTH_START_RHO; rho 0 is the control's join)
rem   distance / shorten belt nodes that margin shortening thickens keep their lithosphere and thicken it by pure shear
rem                      (bounded at the frame's top); belt nodes not shortened are reset by the arc as before
rem   crust / shorten    both
rem  Rule. Vetoes, every one on its arm: V1 pytest passes (incl. test_lt_link); V2' elevation gates +-300 m on 6 of 6
rem   seeds; V3 ledgers close; V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the
rem   highest node over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: both options if crust / shorten passes every veto. Otherwise each option whose single arm passes every
rem   veto; if both single arms pass and the combined arm fails, 'shorten' alone (the larger measured misfit, fixed
rem   here). An option is more realistic by measurement and mechanism whatever a score does; a worse score becomes the
rem   next thing to work on.
rem  Cloud pre-check, seed 1012 (unused), 250 Myr, control and crust / shorten, read only, not pooled (journal 105):
rem   bounded (L <= 177.7 km, highest node 5,976 against 6,214 m); land p50 1,133 -> 1,046 m, p90 2,270 -> 1,957, sea 84 m
rem   higher, dL/dT 0.65 -> 1.21 (the start join); young lithosphere under crust > 45 km stays ~96 km: shortening
rem   thickens 9% of belt nodes per step, the arc reset still fires on 91% (belt_diag.py). Rule unchanged.
rem  Expected (written before the runs, not a band): 'crust' moves the start's dL/dT 0.9-1.3 -> 1.6-1.9 (n = 217
rem   static, seeds 1001 / 1012), little at 0 Myr (land p50 0 to -40 m); 'shorten' raises L under young thick crust
rem   from 97 km toward 120-170, lowers land on crust > 45 km by 100-400 m, land p90 -50 to -150 m, land p50 -10 to
rem   -60 m, dL/dT up.
rem  Scored, reported (no veto): lith_crust2, column_split, where_misfit and young_thick_planet per arm; land fraction
rem   and p25 / p50 / p90, sea, C1, crust bands, export, volume, K1-K4, edge, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P81=%CAL%\precheck-p81-ltlink
%PY% -u -m scripts.run_sweep --name lt-link --grid start_lithosphere=distance,crust --grid belt_lithosphere=reset,shorten --seeds 1001,1002,1003,1004,1005,1006 --control start_lithosphere-distance_belt_lithosphere-reset --identity %CAL%\sweep-20261109-222925\lith_frame-lamb || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P80%\lith_crust2.py "%SW%\start_lithosphere-distance_belt_lithosphere-reset" "%SW%\start_lithosphere-crust_belt_lithosphere-reset" "%SW%\start_lithosphere-distance_belt_lithosphere-shorten" "%SW%\start_lithosphere-crust_belt_lithosphere-shorten" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
for %%a in (distance_belt_lithosphere-reset crust_belt_lithosphere-reset distance_belt_lithosphere-shorten crust_belt_lithosphere-shorten) do (
  %PY% %P75%\column_split.py "%SW%\start_lithosphere-%%a" > "%SW%\reviews\column_split_%%a.txt" 2>&1 || goto fail
  %PY% %P81%\where_misfit.py "%SW%\start_lithosphere-%%a" > "%SW%\reviews\where_misfit_%%a.txt" 2>&1 || goto fail
  %PY% %P81%\where_misfit.py "%SW%\start_lithosphere-%%a" t-0250.0 >> "%SW%\reviews\where_misfit_%%a.txt" 2>&1 || goto fail
  %PY% %P81%\young_thick_planet.py "%SW%\start_lithosphere-%%a" > "%SW%\reviews\young_thick_%%a.txt" 2>&1 || goto fail
)
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
