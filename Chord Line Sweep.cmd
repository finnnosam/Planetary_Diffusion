@echo off
rem Chord Line Sweep (8 Nov; journal 95, seventy-second pass). Rule fixed 8 Nov before the runs.
rem  The choosing reading (the Rift Line Sweep, sweep-20261107-213619, read in pass 72; rift_line 'weak' adopted):
rem  with the split line in weak lithosphere K1 rose 103.9k -> 117.3k km at 0 Myr (window 104.7k -> 113.6k; Earth
rem  173.1k) and the pieces of 1-10 Mkm2 doubled (1 -> 2 per snapshot, 19.1k -> 33.3k km of edge), but pieces under
rem  1 Mkm2 are still absent (planet 2 crumbs, 0.1k km passive; Earth 28 pieces, 17.4k km passive: Madagascar,
rem  Rockall, Seychelles, Lord Howe, Flemish Cap, ...). The cause: every split's line is the bisector of two seeds far
rem  apart in the plate, so it halves the plate (new-plate area p50 11-13 Mkm2) and cuts across the whole continent.
rem  Earth: a rift nucleates at a point and propagates both ways to the plate's free edge (Courtillot 1982, Tectonics 1,
rem  239; Van Wijk and Blackman 2005, EPSL 229, 247); plate forces act at the plate's edges, so the cut through the
rem  nucleus that needs the least force breaks first; most microcontinents are rifted off margins, many beside
rem  subduction (van den Broek and Gaina 2020, Tectonics 39).
rem  This sweep: rift_line 'chord': the nucleus is a random continental node of the plate; the cut is the cheapest pair
rem  of cheapest paths from it to the plate's edge (cost as 'weak': length x wiggle / strength) whose ends lie more than
rem  90 degrees apart as seen from it and that leaves at least piece_min_area_mkm2 (1 Mkm2, the model's smallest plate)
rem  on both sides; if none qualifies, the 'weak' bisector. Ocean-only splits keep the bisector. No new number (5-degree
rem  direction bins, a discretisation). Same rate law (K: Earth's margin per Mkm2, per km of the candidate's margin,
rem  so shorter cuts split more often and the margin made per Mkm2 is kept), same opening (3.5 cm/yr).
rem  Expected (written before the runs, not a band): smaller new plates and continental pieces: pieces under 1 Mkm2
rem  and 1-10 Mkm2 up, the edge per Mkm2 and K1 up; more splits and plates. Risks: slivers thinned below 20 km by the
rem  100 km rift zones (lost to the label-free count), more young floor (T4 up), more microplates merging back.
rem  Pre-check (cloud seed 1007, 250 Myr, read only, not pooled; precheck-p72-chord\seed1007): see journal 95.
rem  Arms (2 x 6 seeds = 12 runs, about 40-50 min); folders rift_line-...:
rem   weak   today's PARAMS (control; must equal rift_line-weak of sweep-20261107-213619: identity.txt)
rem   chord  + the least-force cut through a continental nucleus
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers
rem   close (tests); V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the highest node
rem   over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: rift_line 'chord' if its arm passes every veto. Its case is Earth's propagating rift and the least-force
rem   cut, replacing an assumed line that halves its plate; not a score.
rem  Scored, reported (no veto): the edge by piece size (piece_edge: Earth 29 pieces under 10 Mkm2, 26k km passive),
rem   K1 at 0 Myr (edge_split) and over the window (edge_window), the unbroken-rift edge (ridge_near2; Earth CRB 14.8k)
rem   and its rift pairs (rift_pairs), the chord share and new-plate areas (chord_events), K2-K4 (margin_age, w_read),
rem   margins band, rifted production, the rift and survival records, land p25 / p50 / p90, coast_budget, export,
rem   volume, label-free area, T2-T5, plates.
rem  Runs with scripts.run_sweep: tests, the 12 runs and the reviews under one progress line; then the edge probes
rem  below, each into <sweep>\reviews\*.txt.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P67=%CAL%\precheck-p67-survival
set P69=%CAL%\precheck-p69-gate
set P70=%CAL%\precheck-p70-edge
set P71=%CAL%\precheck-p71-window
set P72=%CAL%\precheck-p72-chord
%PY% -u -m scripts.run_sweep --name chord-line --grid rift_line=weak,chord --seeds 1001,1002,1003,1004,1005,1006 --control rift_line-weak --identity %CAL%\sweep-20261107-213619\rift_line-weak || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
set A=%SW%\rift_line-weak
set B=%SW%\rift_line-chord
echo Edge probes on %SW%
%PY% %P69%\edge_split.py "%A%" "%B%" > "%SW%\reviews\edge_split.txt" 2>&1 || goto fail
%PY% %P71%\edge_window.py "%A%" "%B%" > "%SW%\reviews\edge_window.txt" 2>&1 || goto fail
%PY% %P71%\piece_edge.py "%A%" "%B%" > "%SW%\reviews\piece_edge.txt" 2>&1 || goto fail
%PY% %P72%\chord_events.py "%A%" "%B%" > "%SW%\reviews\chord_events.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%A%" > "%SW%\reviews\ridge_near2_weak.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%B%" > "%SW%\reviews\ridge_near2_chord.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%A%" > "%SW%\reviews\rift_pairs_weak.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%B%" > "%SW%\reviews\rift_pairs_chord.txt" 2>&1 || goto fail
%PY% %P67%\margin_age_noearth.py "%A%\seed1001" "%A%\seed1002" "%A%\seed1003" "%A%\seed1004" "%A%\seed1005" "%A%\seed1006" > "%SW%\reviews\margin_age_weak.txt" 2>&1 || goto fail
%PY% %P67%\margin_age_noearth.py "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\margin_age_chord.txt" 2>&1 || goto fail
%PY% %P67%\w_read.py "%A%\seed1001" "%A%\seed1002" "%A%\seed1003" "%A%\seed1004" "%A%\seed1005" "%A%\seed1006" "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\w_read.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
