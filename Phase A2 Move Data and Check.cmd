@echo off
rem Phase A2 of the migration plan (steps 5-6): move data and the planet chain's outputs here from
rem ..\Infinite_Diffusion (a rename on the same drive), run the tests and the identity check, then push.
rem Run no sweep in either folder while this runs. The identity run takes about 25 minutes.
cd /d "%~dp0"
set "OLD=%~dp0..\Infinite_Diffusion"
set "PY=.venv\Scripts\python.exe"
set "CAL=outputs\planet-calibration"
set "REF=%CAL%\sweep-20261028-200644\trench_polarity-node_trench_flux-band_arc_H_ocean_m-300_arc_accretion-False\seed1001"
if not exist %PY% (
  echo Run Phase A1 Import History.cmd first.
  goto fail
)

echo === 0/4 Commit the import fixes: LF line endings, six missed planet_training modules, pytest temp folder ===
git add .gitattributes pyproject.toml planet_training "Phase A2 Move Data and Check.cmd" || goto fail
git add --renormalize . || goto fail
git diff --cached --quiet || git commit -F migration\phase-a-fix-commit-msg.txt || goto fail

echo === 1/4 Move data and the planet chain's outputs ===
if not exist data (
  move "%OLD%\data" data >nul || goto fail
)
if not exist outputs mkdir outputs
for %%d in (planet-calibration planet-landscape planet-landscape-timing planet-meso-check planet-plates planet-reference erosion-law-check meso meso-plains-diagnosis toy_meso tectonic-prototype claude-sync) do (
  if exist "%OLD%\outputs\%%d" if not exist "outputs\%%d" (
    echo   outputs\%%d
    move "%OLD%\outputs\%%d" "outputs\%%d" >nul || goto fail
  )
)

echo === 2/4 All moved tests ===
%PY% -m pytest -q tests > outputs\phase-a-pytest.txt 2>&1
set RC=%errorlevel%
type outputs\phase-a-pytest.txt
if not "%RC%"=="0" goto fail

echo === 3/4 Identity: seed 1001, today's model, against the Trench Flux sweep's control ===
%PY% -u -m scripts.sweep_plates --rule shelf --seeds 1001 --grid erosion_K_factor,erosion_theta,erosion_Sc,margin_shortening,sediment_accretion,continental_sediment,crust_kappa_normal_m2_yr,sediment_compaction,piece_rule=1.0:200:0.036:True:0.13:True:30:True:plates --grid trench_polarity,trench_flux,arc_H_ocean_m,arc_accretion=node:band:300:False || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
%PY% tools\compare_runs.py "%SW%\combo-01\seed1001" "%REF%" > outputs\phase-a-identity.txt 2>&1
type outputs\phase-a-identity.txt
findstr /b /c:"IDENTICAL" outputs\phase-a-identity.txt >nul || goto fail

echo === 4/4 Push to GitHub ===
git push origin main || goto fail
echo.
echo Phase A2 finished: tests pass, the run is IDENTICAL, and the repository is on GitHub.
pause
exit /b 0
:fail
echo.
echo *** Stopped: the step above failed. Nothing after it ran. Copy the message above to Claude. ***
pause
exit /b 1
