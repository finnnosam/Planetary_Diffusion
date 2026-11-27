@echo off
rem Render Surface (20 Nov; journal 127). A full-resolution elevation map (topography colours, NW hillshade) of a
rem  surface, written beside its source: <seed>\surface-t+0000.0.png, or <landscape run>\surface.png.
rem  Drag onto this file: a seed folder (its 0 Myr snapshot), an arm or sweep folder (each arm's first seed),
rem  a snapshot .npz, or a landscape run folder (outputs\planet-landscape\..., the 11.53 km surface).
rem  Double-click: the newest sweep, each arm's first seed, at 0 Myr.
rem  Native resolution: one pixel per grid node at the equator, nearest-node values, nothing interpolated.
rem  Options for scripts.render_surface (PowerShell): --time -100, --all-nodes (1.5x pixel grid: every node drawn), --boundaries.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set "TARGET=%~1"
if not "%TARGET%"=="" goto go
for /f "delims=" %%d in ('dir /b /ad /o-n outputs\planet-calibration\sweep-*') do (set "TARGET=outputs\planet-calibration\%%d" & goto go)
:go
echo Rendering %TARGET%
%PY% -m scripts.render_surface "%TARGET%" || goto fail
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
