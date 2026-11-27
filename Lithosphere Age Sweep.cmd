@echo off
rem Lithosphere Age Sweep (1 Nov; journal 78, fifty-sixth pass). Item 1 (continental crust distribution): C, the
rem  continental lithosphere's thermal state (planet_tectonics.lithosphere; PARAMS lithosphere 'thermal_age').
rem  Pass 55 (journal 77): with the collapse flow off the thick crust takes Earth's shape, but the continent is still a
rem  dome (P5 fails everywhere) and export fails. From coast to interior Earth's land falls 773 m below the planet's
rem  Airy line (old, cold mantle lithosphere); the planet stands on it everywhere (+206 m). C carries a tectonothermal
rem  age per continental node (start: Earth's CRUST1.0 types joined by the rank of distance from the continent's edge;
rem  reset by collision contacts, continental arcs over converging trenches and new continental crust; thinned with
rem  the crust by stretching) and
rem  stands the surface at -(30 / 3300)(L - 130.67 km), L = 93.6 + 0.04 t km (Artemieva 2006; Lamb et al. 2020).
rem  No fitted number.
rem  Arms (3 x 6 seeds = 18 runs, about 1.5-2 h), everything else PARAMS; folders combo-01..03 (params.json):
rem   combo-01 crust_kappa_m2_yr 0, trench_polarity pair, trench_flux consumed, arc_H_ocean_m 1000, arc_accretion True,
rem       lithosphere none: combo-03 of sweep-20261101-115244 rerun; must equal it field by field (identity.txt)
rem   combo-02 crust_kappa_m2_yr 0, the rest as today's PARAMS, lithosphere thermal_age (combo-02 + C)
rem   combo-03 combo-01's settings with lithosphere thermal_age (combo-03 + C)
rem  Rule (fixed 1 Nov before the runs), six-seed means unless stated:
rem   C1 coast-to-interior fall of z - Airy (coast_budget: 600-1200 km minus 0-50 km inland) in 0.5-1.5 x Earth's
rem      (Earth -773 m: -387 ... -1,160 m)
rem   C2 = P5 rise of land inland (coast_profile: 600-1200 minus 0-50 km) <= 2 x Earth's (Earth 386 m: 772 m)
rem   P1 > 45 km share of the label-free continent <= 0.13; P2 > 55 km <= 0.036; P3 slope of land > 2 km p50 in
rem      7.0-17.6 m/km; P4 > 45 km depth from its edge p50 in 58-146 km (crust_review, interior_probe; as the
rem      Collapse Flow Sweep)
rem   Guards (all must hold): pytest passes; elevation gates +-300 m on 6 of 6 seeds (elevation_check); continental
rem   volume change >= -14% (crust_budget); export 4-8 km3/yr (crust_budget); margins pass (margin_width); plates at
rem   0 Myr <= 40 on every seed; ledgers close (tests).
rem   An arm is adopted if C1 and C2 pass, >= 3 of P1-P4 pass and every guard holds. If both C arms qualify, combo-03
rem   (it carries the volume source). Adopting sets PARAMS lithosphere 'thermal_age' and crust_kappa_m2_yr 0 (and, for
rem   combo-03, B-d2 arm C's settings). If neither: nothing adopted; the failing test or guard sets the next step:
rem   C1 fails but C1 without the planet's dynamic topography (coast_budget, reported) is in the band -> the dynamic
rem   topography from Earth's measured spectrum next (STATE next item 4), then C is read again; C1 fails both ways ->
rem   where the term's fall is lost (ages and resets by band); export fails -> the erosion rate of narrow belts.
rem  Pre-check (cloud seed 1007, combo-03 + C against combo-03, journal 78): C1 +247 -> +112 m (fails); without the
rem  dynamic topography +39 -> -233 m. The random field (fitted) puts interior land on its highs (+150-240 m).
rem  Reported, not scored: land p25 / p50 / p90, land fraction, continental crust drowned, area, sea level on the datum,
rem  lowland beyond 300 km, lithosphere age by band and C1 without the dynamic topography (coast_budget), T2-T5
rem  (plate_pieces_review), arc growth A1.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py tests\test_lithosphere.py > %CAL%\lithosphere-age-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\lithosphere-age-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid crust_kappa_m2_yr,trench_polarity,trench_flux,arc_H_ocean_m,arc_accretion,lithosphere=0:pair:consumed:1000:True:none,0:node:band:300:False:thermal_age,0:pair:consumed:1000:True:thermal_age || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\lithosphere-age-pytest.txt "%SW%\reviews\pytest.txt" >nul
%PY% tools\compare_runs.py "%SW%\combo-01" %CAL%\sweep-20261101-115244\combo-03 > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% -u -m planet_reviews.coast_budget "%SW%" %CAL%\sweep-20261101-115244\combo-02 --json "%SW%\reviews\coast_budget.json" > "%SW%\reviews\coast_budget.txt" 2>&1
type "%SW%\reviews\coast_budget.txt"
%PY% -u -m planet_reviews.coast_profile "%SW%" --json "%SW%\reviews\coast_profile.json" > "%SW%\reviews\coast_profile.txt" 2>&1
type "%SW%\reviews\coast_profile.txt"
%PY% -u -m scripts.plate_pieces_review --latest --control combo-01 > "%SW%\reviews\plate_pieces_review.txt" 2>&1
type "%SW%\reviews\plate_pieces_review.txt"
for %%r in (crust_review crust_budget elevation_check margin_width interior_probe) do (
  %PY% -u -m scripts.%%r --latest > "%SW%\reviews\%%r.txt" 2>&1
  type "%SW%\reviews\%%r.txt"
)
echo Reviews saved in %SW%\reviews
echo done
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
