@echo off
rem Initiation Site Sweep (6 Nov; journal 89, sixty-seventh pass). Rule fixed 6 Nov before the runs.
rem  Read of the Start Breakup Sweep (sweep-20261106-211658): V1-V5 pass on every arm (pytest 135; gates 6 of 6; T1
rem  <= 0.0042; highest node <= 8,195 m); S1 (rifted production -240 ... -200 Myr <= 1,000 km/Myr): whole/coast/forces
rem  198, whole/voronoi 638 pass, whole/coast/random 1,434 and the control 2,553 fail. By the rule's order combo-04
rem  (start_breakup whole, start_plates coast, start_poles forces) is adopted: it is PARAMS now (the bare run is
rem  identical step by step to its settings on the old code, n = 48, 2 seeds x 70 steps, threads 1). Scored: volume
rem  -0.8%, label-free area 0.357, crust median 38.9 km, > 45 km 0.169, interior crust 40 km (Earth 39.5-41), land p50
rem  1,295 m, export 1.4 km3/yr (4-8), K1 96.2k km, K2 62.3k, K3 p50 75 Myr (fails, 85-141), K4 15.2 Mkm2 / 59 km2.
rem  Survival record, last 150 Myr: lifespan 174 Myr (Earth 130-181), births 428 km/Myr (Earth's steady state needs
rem  ~960-1,330): the passive edge is short by births. Of its continental edge with margin age >= 50 Myr, 0.53 is
rem  within 300 km of a boundary (Earth 0.21); the nearest a trench at 0.56 of it.
rem  Built this pass (default off / output only; the new PARAMS path identical step by step):
rem   init_site ocean: spontaneous initiation on a plate that carries continent no longer breaks its ocean off along
rem     the coast; old floor welded to the continent does not initiate, old floor apart from it breaks off where it
rem     meets younger floor. Earth: no unambiguous passive margin turned convergent in the Cenozoic (Stern 2004; Stern &
rem     Gerya 2018); old margins too strong to fail unforced (Mueller & Phillips 1991). No new number. Skips are
rem     recorded as events 'initiation at margin skipped'. On the adopted arm 1 of 15 initiations split a plate with
rem     continent: expected small.
rem   planet_reviews.passive_survival.RiftRecorder (stats 'rift_record', every run): young rifted margin (flagged edge beside a
rem     ridge and floor < 10 Myr) carried with its plate until it becomes passive edge or ends (subduction, collision,
rem     closed, lost, other, stalled after 100 Myr). Earth's steady state: ~0.6-0.85 of rifted margin becomes passive.
rem  Arms (4 x 6 seeds = 24 runs, about 60-75 min); folders init_site-..._subduction_initiation-...:
rem   margin:instant      today's PARAMS (control; must equal start_breakup-whole_start_plates-coast_start_poles-forces
rem                       of sweep-20261106-211658 on the shared values: identity.txt)
rem   ocean:instant       + init_site ocean
rem   margin:convergence  + subduction_initiation convergence (B-d3, built 30 Oct, never swept: a slab only after
rem                       slab_init_km 125 of underthrusting, Gurnis et al. 2004; Bird's OCB 17.4k km of slow
rem                       ocean-ocean convergence without a slab)
rem   ocean:convergence   both
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2 elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close
rem   (tests); V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the highest node over the
rem   run <= 8,850 m or no higher than the control's there.
rem  B-d3's own test O1: ocean-ocean convergence without a slab (plate_pieces_review, window -100 ... 0 Myr, six-seed
rem   mean) >= 8.7k km (half of Bird's OCB 17.4k).
rem  Adoption, two decisions on the 2 x 2: D1 init_site ocean if ocean:instant passes every veto (its case is Earth's
rem   record, not a score); D2 subduction_initiation convergence if margin:convergence passes every veto and O1. If both,
rem   ocean:convergence must pass every veto and O1 too, else only D1. Adopted settings become PARAMS.
rem  Scored, reported (no veto): the rift record (share of young rifted margin becoming passive, its fates, time to
rem   passive), survival (births, lifespan, fates), skipped initiations, T2-T5, and read after by Claude: the birth
rem   probe (passive share of 50-200 Myr margins, Earth 0.79), K1-K4; export, volume, margins band, land p25 / p50 /
rem   p90, > 45 / > 55 km, label-free area, sea on the datum, C1, foundered volume.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py tests\test_lithosphere.py tests\test_dyn_measured.py tests\test_margin_wedge.py tests\test_collision_spread.py tests\test_breakup_match.py tests\test_start_breakup.py tests\test_init_site.py > %CAL%\init-site-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\init-site-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid init_site,subduction_initiation=margin:instant,ocean:instant,margin:convergence,ocean:convergence || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
set CTL=%SW%\init_site-margin_subduction_initiation-instant
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\init-site-pytest.txt "%SW%\reviews\pytest.txt" >nul
%PY% tools\compare_runs.py "%CTL%" %CAL%\sweep-20261106-211658\start_breakup-whole_start_plates-coast_start_poles-forces > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% tools\max_elevation.py "%SW%" > "%SW%\reviews\max_elevation.txt" 2>&1
type "%SW%\reviews\max_elevation.txt"
%PY% -u -m planet_reviews.rift_production "%SW%" --json "%SW%\reviews\rift_production.json" > "%SW%\reviews\rift_production.txt" 2>&1
type "%SW%\reviews\rift_production.txt"
%PY% -u -m planet_reviews.passive_survival "%SW%" --json "%SW%\reviews\passive_survival.json" > "%SW%\reviews\passive_survival.txt" 2>&1
type "%SW%\reviews\passive_survival.txt"
%PY% -u -m planet_reviews.coast_budget "%SW%" "%CTL%" --json "%SW%\reviews\coast_budget.json" > "%SW%\reviews\coast_budget.txt" 2>&1
type "%SW%\reviews\coast_budget.txt"
%PY% -u -m planet_reviews.coast_profile "%SW%" "%CTL%" --json "%SW%\reviews\coast_profile.json" > "%SW%\reviews\coast_profile.txt" 2>&1
type "%SW%\reviews\coast_profile.txt"
%PY% -u -m scripts.plate_pieces_review --latest --control init_site-margin_subduction_initiation-instant > "%SW%\reviews\plate_pieces_review.txt" 2>&1
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
