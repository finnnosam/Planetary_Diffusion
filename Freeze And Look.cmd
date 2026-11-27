@echo off
rem Freeze and look (31 Oct; journal 74, re-evaluation step 3, agreed with the user). The whole chain once on today's
rem  baseline (PARAMS, no overrides), seed 1001: plates (250 Myr, saved every Myr over the last 40, for the landscape
rem  stage) -> landscape (11.53 km, last 30 Myr, theta 5) -> 23 km hand-off and meso against Earth's windows (GPU).
rem  Nothing is adopted from it: it is read for what the terrain gets wrong now (hypsometry against Earth, maps at
rem  65,536 colour levels, the meso's flat-window check), to set the order of the next items.
rem  The last landscape and meso runs were on the 22 Oct baseline (outputs\planet-landscape\theta-20261022-*).
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set PL=outputs\planet-plates\freeze-20261031
set LS=outputs\planet-landscape\freeze-20261031-seed1001
set MC=outputs\planet-meso-check\freeze-20261031-seed1001
%PY% -m pytest -q tests\test_shelf.py tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py || goto fail
%PY% -u -m scripts.simulate_plates --seeds 1001 --out %PL% || goto fail
%PY% -u -m scripts.run_planet_landscape --plates %PL%\seed1001 --out %LS% || goto fail
%PY% -u -m scripts.planet_meso_check --landscape %LS% --out %MC% || goto fail
echo Plates: %PL%  Landscape: %LS%  Meso check: %MC%
echo done
pause
exit /b 0
:fail
echo FAILED
pause
exit /b 1
