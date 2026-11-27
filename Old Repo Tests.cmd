@echo off
rem Phase C follow-up, second run (local working script, not published). The first run staged the removal of the
rem tests and scripts for branch-only code; two remaining tests passed --keep-checkpoints, a branch-only flag, and
rem Claude has removed it from them. Run the generator's whole test suite, then commit and push.
cd /d "%~dp0..\Infinite_Diffusion"
set "NEW=%~dp0"
set "OPY=.venv\Scripts\python.exe"

echo === 1/3 The changes to commit ===
git add tests/test_cube_workflow.py tests/test_draft.py || goto fail
git diff --cached --stat
git status --porcelain | findstr /v /b /c:"D  " /c:"M  " > "%NEW%outputs\phase-c3-other.txt"
for %%A in ("%NEW%outputs\phase-c3-other.txt") do if %%~zA gtr 0 (
  type "%NEW%outputs\phase-c3-other.txt"
  echo Other uncommitted changes in Infinite_Diffusion; nothing was committed.
  goto fail
)

echo === 2/3 The generator's whole test suite ===
%OPY% -m pytest -q -p no:cacheprovider --basetemp=outputs\pytest-tmp-phase-c3 tests > "%NEW%outputs\phase-c3-pytest.txt" 2>&1
set RC=%errorlevel%
type "%NEW%outputs\phase-c3-pytest.txt"
if not "%RC%"=="0" goto fail

echo === 3/3 Commit and push ===
git commit -q -F "%NEW%migration\phase-c2-commit-msg-v2.txt" || goto fail
git push origin main || goto fail
echo.
echo Done: the generator's tests pass on main, and it is pushed.
pause
exit /b 0

:fail
echo.
echo *** Stopped: the step above failed. Nothing was committed. Copy the message above to Claude. ***
pause
exit /b 1
