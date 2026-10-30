@echo off
rem Phase A1 of the migration plan (steps 3-4): import the planet chain from ..\Infinite_Diffusion with its
rem history, and create this repository's environment. Needs Phase 0 done (tag pre-split-2026-10-29).
rem Pushes nothing; Phase A2 pushes after the identity check. Double-click from this folder.
cd /d "%~dp0"
set "OLD=%~dp0..\Infinite_Diffusion"
set "SPLIT=%TEMP%\pd-split-import"
set "PY=.venv\Scripts\python.exe"
git -C "%OLD%" rev-parse -q --verify refs/tags/pre-split-2026-10-29 >nul
if errorlevel 1 (
  echo Phase 0 has not run yet: the tag pre-split-2026-10-29 is missing in Infinite_Diffusion.
  goto fail
)

echo === 1/5 Environment: same Python and package versions as Infinite_Diffusion (torch is a large download) ===
if exist .venv\Scripts\python.exe goto haveenv
"%OLD%\.venv\Scripts\python.exe" -m venv .venv || goto fail
:haveenv
"%OLD%\.venv\Scripts\python.exe" -m pip freeze --exclude-editable > "%TEMP%\pd-requirements.txt" || goto fail
%PY% -m pip install --upgrade pip || goto fail
%PY% -m pip install -r "%TEMP%\pd-requirements.txt" --extra-index-url https://download.pytorch.org/whl/cu130 || goto fail
%PY% -m pip install h5py git-filter-repo || goto fail
%PY% -m pip install -e "%OLD%" --no-deps || goto fail

echo === 2/5 Filtered copy of the old history (only the moving files) ===
if exist "%SPLIT%" rmdir /s /q "%SPLIT%"
git clone --no-local "%OLD%" "%SPLIT%" || goto fail
git -C "%SPLIT%" reset --hard pre-split-2026-10-29 || goto fail
pushd "%SPLIT%"
"%~dp0.venv\Scripts\git-filter-repo.exe" --paths-from-file "%~dp0migration\split-paths.txt" --force
if errorlevel 1 (
  popd
  goto fail
)
popd

echo === 3/5 Merge it into this repository ===
git fetch --no-tags "%SPLIT%" main:split-import || goto fail
git merge --allow-unrelated-histories --no-edit -m "Import the planet chain from Infinite_Diffusion with its history (tag pre-split-2026-10-29)" split-import || goto fail
git branch -D split-import

echo === 4/5 Repository files ===
git add .gitignore pyproject.toml README.md migration "Phase A1 Import History.cmd" "Phase A2 Move Data and Check.cmd" || goto fail
git commit -F migration\phase-a-commit-msg.txt || goto fail

echo === 5/5 Quick tests (no data needed) ===
%PY% -m pytest -q tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_shelf.py
if errorlevel 1 goto fail
echo.
git log --oneline -3
echo Phase A1 finished. Next: Phase A2 Move Data and Check.cmd (run no sweep while it runs).
pause
exit /b 0
:fail
echo.
echo *** Stopped: the step above failed. Nothing after it ran. Copy the message above to Claude. ***
pause
exit /b 1
