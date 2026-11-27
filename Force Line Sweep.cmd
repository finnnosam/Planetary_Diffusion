@echo off
rem Force Line Sweep (8 Nov; journal 96, seventy-third pass). Rule fixed 8 Nov before the runs.
rem  The choosing reading (the Chord Line Sweep, sweep-20261108-092844, read in pass 73; rift_line 'chord' adopted by
rem  its rule): the chord makes Earth's small pieces (0.1-1 Mkm2: 0 -> 2 per snapshot, 10.5k km of edge; Earth 3 pieces,
rem  9.2k km), new plates p50 12.9 -> 2.1 Mkm2, plates at 0 Myr 29 -> 44 (Earth 52), rifted production 1.30x -> 1.11x
rem  Earth, land p50 1,228 -> 1,083 m; but the passive edge halves on the big pieces (99.6k -> 51.9k km; K1 at 0 Myr
rem  117.3k -> 74.1k, window 113.6k -> 67.2k; Earth 173.1k), lifespan 173 -> 126 Myr, T2 28.8k -> 40.0k, T4 4.40 -> 4.85.
rem  The cause: the cut that breaks is the weakest through a random nucleus, and the weakest is a short cut across a
rem  corner: the big plate (plate 0 in most rift pairs) sheds slivers instead of splitting. On Earth a cut breaks when
rem  the force it must transmit exceeds its strength (Forsyth and Uyeda 1975 for the plate forces; Buck 1991 and
rem  Kusznir and Park 1987 for rifting when tension exceeds the lithosphere's integrated strength): a corner cut
rem  carries almost no force.
rem  This sweep: rift_line 'force': the same nucleus and candidate cuts as 'chord'; the one that breaks is the one
rem  nearest failure, its opening force over its cost (the cut's strength). The force: the plate moving as one under its
rem  boundary forces and drag as force_law 'torque' books them (slab pull 5e13, ridge push 3e12 on oceanic ridges, slab
rem  suction 0, basal drag with continental keels 3.4x, slab resistance at 8.4 cm/yr: its fit to Earth's plates), the
rem  torque the cut carries about the pole that opens it (split()'s). Candidates ranked on the plate's sectors as seen
rem  from the nucleus, checked exactly in that order. No candidate in tension (a plate nothing pulls apart: an
rem  overriding continent under the fit's zero suction, e.g. the starting supercontinent): the cheapest cut, as 'chord'.
rem  No new number. Same rate law (K) and opening (3.5 cm/yr). Collisions' coupling to neighbours left out (assumed).
rem  Expected (written before the runs, not a band): cuts on plates with slabs larger and across them (fewer slivers,
rem  K1 and the big pieces' passive edge up, toward the weak arm's), the supercontinent unchanged (no force on it);
rem  maybe fewer small pieces. Pre-check (cloud seed 1007, read only, not pooled; precheck-p73-force\seed1007): journal 96.
rem  Arms (2 x 6 seeds = 12 runs, about 40-55 min); folders rift_line-...:
rem   chord  today's PARAMS (control; must equal rift_line-chord of sweep-20261108-092844: identity.txt)
rem   force  + the cut nearest failure (opening force over strength)
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers
rem   close (tests); V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the highest node
rem   over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: rift_line 'force' if its arm passes every veto. Its case is the failure criterion (a cut breaks where the
rem   force it carries is largest for its strength), replacing a criterion of strength alone; not a score.
rem  Scored, reported (no veto): the edge by piece size (piece_edge), K1 at 0 Myr (edge_split) and over the window
rem   (edge_window), the unbroken-rift edge and rift pairs, force / chord / bisector shares and new-plate areas
rem   (chord_events), K2-K4 (margin_age, w_read), margins band, rifted production, the rift and survival records, land
rem   p25 / p50 / p90, coast_budget, export, volume, label-free area, T2-T5, plates.
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
%PY% -u -m scripts.run_sweep --name force-line --grid rift_line=chord,force --seeds 1001,1002,1003,1004,1005,1006 --control rift_line-chord --identity %CAL%\sweep-20261108-092844\rift_line-chord || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
set A=%SW%\rift_line-chord
set B=%SW%\rift_line-force
echo Edge probes on %SW%
%PY% %P69%\edge_split.py "%A%" "%B%" > "%SW%\reviews\edge_split.txt" 2>&1 || goto fail
%PY% %P71%\edge_window.py "%A%" "%B%" > "%SW%\reviews\edge_window.txt" 2>&1 || goto fail
%PY% %P71%\piece_edge.py "%A%" "%B%" > "%SW%\reviews\piece_edge.txt" 2>&1 || goto fail
%PY% %P72%\chord_events.py "%A%" "%B%" > "%SW%\reviews\chord_events.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%A%" > "%SW%\reviews\ridge_near2_chord.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%B%" > "%SW%\reviews\ridge_near2_force.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%A%" > "%SW%\reviews\rift_pairs_chord.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%B%" > "%SW%\reviews\rift_pairs_force.txt" 2>&1 || goto fail
%PY% %P67%\margin_age_noearth.py "%A%\seed1001" "%A%\seed1002" "%A%\seed1003" "%A%\seed1004" "%A%\seed1005" "%A%\seed1006" > "%SW%\reviews\margin_age_chord.txt" 2>&1 || goto fail
%PY% %P67%\margin_age_noearth.py "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\margin_age_force.txt" 2>&1 || goto fail
%PY% %P67%\w_read.py "%A%\seed1001" "%A%\seed1002" "%A%\seed1003" "%A%\seed1004" "%A%\seed1005" "%A%\seed1006" "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\w_read.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
