@echo off
rem Rift Line Sweep (7 Nov; journal 94, seventy-first pass). Rule fixed 7 Nov before the runs.
rem  The choosing reading (precheck-p71-window, on the Rift Speed Sweep sweep-20261107-191259, window means over the
rem  saved snapshots -40 ... 0 Myr): with the fast rift phase (adopted, pass 71) the unbroken-rift edge is 19.5k km
rem  (Earth CRB-near 14.8k; the 0 Myr snapshot read 30.3k: three seeds caught fresh long splits) and the non-passive
rem  edge is Earth's length (~139k against 131k km). K1's gap (passive 105k against 173k) is now the edge's length:
rem  Earth's label-free continent (CRUST1.0 1 degree, n = 217) is one 176 Mkm2 piece + Antarctica (17 Mkm2, 22.7k km of
rem  passive edge) + 29 pieces under 10 Mkm2, mostly rifted microcontinents (Zealandia, Madagascar, Rockall, Seychelles,
rem  Lord Howe, Flemish Cap; oceanic plateaus a minor share) with 26k km of passive edge; the planet's pieces under
rem  10 Mkm2 carry 8k, and its big pieces 1,274 km of edge per Mkm2 against Earth's 1,412. The planet's split line is
rem  the cost-distance bisector of two random seeds far apart in the plate: every split halves its plate, and the belt
rem  term lowers the cost on belts, which moves the line off them (the code comment said it "prefers old sutures").
rem  Earth: continents deform and break where their lithosphere is weak (Audet and Buergmann 2011, Nature Geosci. 4,
rem  184: elastic thickness > 100 km in cratons, < 40 km in Phanerozoic belts and along most margins; Buiter and Torsvik
rem  2014: Atlantic rifts along old belts).
rem  This sweep: rift_line 'weak', the split line's cost per edge is its length / strength, strength = C's thermal
rem  lithosphere thickness L(t) / Earth's continental mean (lith_ref_km) on continental nodes (1 elsewhere) x the as-built
rem  belt factor (1 - suture_weakness on belts): weak lithosphere is slow for both fronts, so they meet in it. Strength in
rem  proportion to L is assumed; no fitted number. Same seeds, same rate (K), same opening (rift_speed_cm_yr 3.5).
rem  Expected (written before the runs, not a band): lines in young belts and the young coast lithosphere: more pieces
rem  under 10 Mkm2 and a longer edge per Mkm2, K1 up; the risk is slivers too thin to live (crumbs welded back).
rem  Arms (2 x 6 seeds = 12 runs, about 35-45 min); folders rift_line-...:
rem   voronoi  today's PARAMS (control; must equal rift_speed_cm_yr-3.5 of sweep-20261107-191259: identity.txt)
rem   weak     + the split line in weak lithosphere
rem  Pre-check (cloud seed 1007, 250 Myr, read only, not pooled; precheck-p71-window\seed1007): unbroken-rift edge
rem   over -40..0 Myr 39.6k -> 29.8k km (0 Myr 38.2k -> 19.5k); passive edge over the window 89.5k -> 91.7k (0 Myr
rem   102.3k -> 111.6k); continental edge over the window 282.7k -> 233.0k; pieces under 10 Mkm2 (window) 9 crumbs + 2
rem   pieces of 1-10 Mkm2 (5.9k km passive) -> one piece of 1-10 (7.9k): no fragments on this seed (merges 2 -> 9);
rem   highest node 6,663 / 8,014 m; volume -0.4 / -1.5%; land p50 1,236 / 1,349 m. Against the expectation on one seed
rem   (the trajectories part: 35 / 28 plates over 1 Mkm2); the six decide.
rem  Rule. Vetoes, every one on its arm: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers
rem   close (tests); V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the highest node
rem   over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: rift_line 'weak' if its arm passes every veto. Its case is Earth's measured weakness control of where
rem   continents break, replacing a line that halves its plate and steers off belts; not a score.
rem  Scored, reported (no veto): K1 and the edge split at 0 Myr (edge_split) and over the window (edge_window), the
rem   edge by piece size (piece_edge: Earth 29 pieces under 10 Mkm2, 26k km passive), the unbroken-rift edge
rem   (ridge_near2; Earth CRB 14.8k) and its rift pairs (rift_pairs), K2-K4 (margin_age, w_read), margins band, rifted
rem   production, the rift and survival records, land p25 / p50 / p90, coast_budget, export, volume, label-free area,
rem   T2-T5, plates.
rem  Runs with scripts.run_sweep: tests, the 12 runs and the reviews under one progress line; then the edge probes
rem  below, each into <sweep>\reviews\*.txt.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P67=%CAL%\precheck-p67-survival
set P69=%CAL%\precheck-p69-gate
set P70=%CAL%\precheck-p70-edge
set P71=%CAL%\precheck-p71-window
%PY% -u -m scripts.run_sweep --name rift-line --grid rift_line=voronoi,weak --seeds 1001,1002,1003,1004,1005,1006 --control rift_line-voronoi --identity %CAL%\sweep-20261107-191259\rift_speed_cm_yr-3.5 || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
set A=%SW%\rift_line-voronoi
set B=%SW%\rift_line-weak
echo Edge probes on %SW%
%PY% %P69%\edge_split.py "%A%" "%B%" > "%SW%\reviews\edge_split.txt" 2>&1 || goto fail
%PY% %P71%\edge_window.py "%A%" "%B%" > "%SW%\reviews\edge_window.txt" 2>&1 || goto fail
%PY% %P71%\piece_edge.py "%A%" "%B%" > "%SW%\reviews\piece_edge.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%A%" > "%SW%\reviews\ridge_near2_voronoi.txt" 2>&1 || goto fail
%PY% %P70%\ridge_near2.py "%B%" > "%SW%\reviews\ridge_near2_weak.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%A%" > "%SW%\reviews\rift_pairs_voronoi.txt" 2>&1 || goto fail
%PY% %P70%\rift_pairs.py "%B%" > "%SW%\reviews\rift_pairs_weak.txt" 2>&1 || goto fail
%PY% %P67%\margin_age_noearth.py "%A%\seed1001" "%A%\seed1002" "%A%\seed1003" "%A%\seed1004" "%A%\seed1005" "%A%\seed1006" > "%SW%\reviews\margin_age_voronoi.txt" 2>&1 || goto fail
%PY% %P67%\margin_age_noearth.py "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\margin_age_weak.txt" 2>&1 || goto fail
%PY% %P67%\w_read.py "%A%\seed1001" "%A%\seed1002" "%A%\seed1003" "%A%\seed1004" "%A%\seed1005" "%A%\seed1006" "%B%\seed1001" "%B%\seed1002" "%B%\seed1003" "%B%\seed1004" "%B%\seed1005" "%B%\seed1006" > "%SW%\reviews\w_read.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
