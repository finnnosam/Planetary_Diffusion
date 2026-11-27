@echo off
rem Mantle Flow Sweep (31 Oct; journal 67, forty-sixth pass). D: dynamic topography from the planet's own slab record.
rem  dyn_rule slabs replaces the random knot field (plate-blind, flat spectrum to degree 30, fitted) by the free-slip
rem  mantle flow of the subducted lithosphere (planet_tectonics.mantle_flow: Earth's radial viscosity sh08 from
rem  data\mantle, the planet's own GDH1 thermal mass, slab dip 58 deg, lower-mantle sinking 12 mm/yr, cut 300 km,
rem  expansivity falling with depth, slab_alpha depth). Design and tests: journal 66; built and pre-checked: journal 67.
rem  Cloud pre-check (seed 1001, one seed only): P-a, P-b, P-c, P-d steps pass; amplitude 623 m (band 110-460) and
rem  flow 4.05 cm/yr (band 0.7-2.7) high, so arm B is expected to fail P-d / P-e unless the seeds differ; it is run
rem  because six seeds decide, and P-h (traction, paired) passed: C(d) +0.20, T2 +0.8k km.
rem  Needs data\mantle\visc.sh08 (already on this machine).
rem  Arms (3 x 6 seeds = 18 runs, about 1.5-2 h), everything else PARAMS:
rem   A dyn_rule-random_force_law-lengths_basal_traction-False : control (today's model); must equal arm A of
rem       sweep-20261028-200644 field by field (identity.txt; stats.json differs only by new keys in params_effective)
rem   B dyn_rule-slabs_force_law-lengths_basal_traction-False  : D
rem   C dyn_rule-slabs_force_law-torque_basal_traction-True    : D with the torque law and traction toward the flow
rem  Rule (fixed 31 Oct before the runs), six-seed means unless stated:
rem   B adopted (dyn_rule slabs) if P-a ... P-e pass (dynamic_topography.txt, last 100 Myr) and P-f holds on 6 of 6
rem   seeds: elevation gates +-300 m (elevation_check), plate count (watch > 40), ledgers. Land fraction, drowned
rem   crust, land p50 and texture are reported, not scored. P-d high: no scale factor (alpha(z) is already in).
rem   P-b fails: debug, do not tune.
rem   C adopted on top (force_law torque, basal_traction True) only if B is adopted, C passes the same, its C(d) of
rem   different plates below 2,000 km (trench_pairs.txt) is >= B's + 0.06, and its T2 (plate_pieces_review) is no
rem   more than B's + 2k km.
rem  Read, not scored (P-g): C(d), rho by class, T2, io trench age (trench_pairs.txt, plate_pieces_review.txt).
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py > %CAL%\mantle-flow-pytest.txt 2>&1
set RC=%errorlevel%
type %CAL%\mantle-flow-pytest.txt
if not "%RC%"=="0" goto fail
%PY% -u -m planet_reviews.mantle_flow_earth > %CAL%\mantle-flow-earth.txt 2>&1
type %CAL%\mantle-flow-earth.txt
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001,1002,1003,1004,1005,1006 --grid dyn_rule,force_law,basal_traction=random:lengths:False,slabs:lengths:False,slabs:torque:True || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
mkdir "%SW%\reviews" 2>nul
copy /y %CAL%\mantle-flow-pytest.txt "%SW%\reviews\pytest.txt" >nul
copy /y %CAL%\mantle-flow-earth.txt "%SW%\reviews\mantle_flow_earth.txt" >nul
set CTRL=dyn_rule-random_force_law-lengths_basal_traction-False
%PY% tools\compare_runs.py "%SW%\%CTRL%" %CAL%\sweep-20261028-200644\trench_polarity-node_trench_flux-band_arc_H_ocean_m-300_arc_accretion-False > "%SW%\reviews\identity.txt" 2>&1
type "%SW%\reviews\identity.txt"
%PY% -u -m planet_reviews.dynamic_topography "%SW%" --json "%SW%\reviews\dynamic_topography.json" > "%SW%\reviews\dynamic_topography.txt" 2>&1
type "%SW%\reviews\dynamic_topography.txt"
%PY% -u -m scripts.plate_pieces_review --latest --control %CTRL% > "%SW%\reviews\plate_pieces_review.txt" 2>&1
type "%SW%\reviews\plate_pieces_review.txt"
%PY% -u -m planet_reviews.trench_pairs "%SW%\%CTRL%" "%SW%\dyn_rule-slabs_force_law-lengths_basal_traction-False" "%SW%\dyn_rule-slabs_force_law-torque_basal_traction-True" --json "%SW%\reviews\trench_pairs.json" > "%SW%\reviews\trench_pairs.txt" 2>&1
for %%r in (crust_budget elevation_check crust_review margin_width) do (
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
