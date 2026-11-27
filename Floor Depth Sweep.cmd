@echo off
rem Floor Depth Sweep (13 Nov; journal 107, eighty-second pass). Rule fixed 13 Nov before the runs.
rem  Why: the land's height is the largest failure (STATE section 0) and half of it is the sea, 494 m below the datum
rem  (journal 106: three_factor.py). The largest part of the low sea is the floor standing deeper at a given age
rem  (-318 m of sea level), and 180 m of that is the law: Earth's floor stands above GDH1 + sediment + crust at every
rem  age beyond 10 Myr. GDH1 (Stein & Stein 1992) was fitted to the North Pacific and the north-west Atlantic.
rem  Measured this pass on the n = 217 graph (isostasy.earth_floor_table): ETOPO on dated floor (Seton 2020, node >= 99%
rem  dated, CRUST1.0 total < 20 km) with the model's own sediment term on GlobSed and crust term on CRUST1.0 taken out,
rem  area means in 26 age bins, pooled monotone: normal oceanic crust (types A1, A0) stands 120-220 m above GDH1
rem  beyond 10 Myr and flattens at 5.35-5.41 km (GDH1 5.65); every dated floor (plateaus in the mean) 130-360 m.
rem  Arms (3 x 6 seeds = 18 runs, about 60 min); folders ocean_depth_rule-...:
rem   gdh1       today's PARAMS (control; must equal start_lithosphere-crust_belt_lithosphere-shorten of
rem              sweep-20261112-203456: identity.txt)
rem   earth      Earth's normal floor by age (the candidate)
rem   earth_all  every dated floor by age: the volume of plateaus and melt-affected crust spread over all floor.
rem              READ ONLY: the planet has no process for thickened oceanic crust and this is no model of one; it
rem              measures what that volume is worth in a run
rem  Rule. Vetoes, every one on its arm: V1 pytest passes (incl. test_ocean_depth); V2' elevation gates +-300 m on 6 of
rem   6 seeds (the gate reads the run's own floor law); V3 ledgers close; V4 plate_pieces_review T1 detached area
rem   <= 0.01 on every seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption: 'earth' if its arm passes every veto (a measurement of Earth's whole normal floor replaces a fit to two
rem   basins; a worse score becomes the next thing to work on). 'earth_all' is not adopted whatever it does.
rem  Cloud pre-check, seed 1013 (unused), control and 'earth', read at -10 Myr (the workspace paused before 0 Myr;
rem   read only, not pooled; journal 107): bounded (highest node 5,863 against 5,882 m, gates -117 / -38 m young
rem   ocean); the sea 128 m higher at -250 Myr, 133 at -100, 105 at -10; land p50 -98 m at -250 Myr.
rem  Expected (written before the runs, not a band; static on today's 0 Myr states, static_floor.py): 'earth': the
rem   sea +124 m on the datum (sd 11), room below the datum 2,973 -> ~2,886 m, land share 0.307 -> ~0.295, land p50
rem   -87 m (sd 17), p90 -99; 'earth_all': the sea +168 m, land p50 -114 m. Between arms on the same seed the L-T Link
rem   Sweep's trajectories differed by 44-126 m of sea and 71-171 m of land p50 (sd), and land p50 followed the sea
rem   0.6-0.7 m per m: read the sea and the room first, land p50 with its seed spread.
rem  Scored, reported (no veto): three_factor (the split), basin (the sea's budget), continent, column_split and
rem   lith_crust2 per arm; land fraction and p25 / p50 / p90, C1, crust bands, export, volume, K1-K4, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
%PY% -u -m scripts.run_sweep --name floor-depth --grid ocean_depth_rule=gdh1,earth,earth_all --seeds 1001,1002,1003,1004,1005,1006 --control ocean_depth_rule-gdh1 --identity %CAL%\sweep-20261112-203456\start_lithosphere-crust_belt_lithosphere-shorten || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P80%\lith_crust2.py "%SW%\ocean_depth_rule-gdh1" "%SW%\ocean_depth_rule-earth" "%SW%\ocean_depth_rule-earth_all" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
for %%a in (gdh1 earth earth_all) do (
  %PY% %P82%\three_factor.py "%SW%\ocean_depth_rule-%%a" > "%SW%\reviews\three_factor_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\basin.py "%SW%\ocean_depth_rule-%%a" > "%SW%\reviews\basin_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\continent.py "%SW%\ocean_depth_rule-%%a" > "%SW%\reviews\continent_%%a.txt" 2>&1 || goto fail
  %PY% %P82%\hyps.py "%SW%\ocean_depth_rule-%%a" > "%SW%\reviews\hyps_%%a.txt" 2>&1 || goto fail
  %PY% %P75%\column_split.py "%SW%\ocean_depth_rule-%%a" > "%SW%\reviews\column_split_%%a.txt" 2>&1 || goto fail
)
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
