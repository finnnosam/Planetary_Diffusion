@echo off
rem Arc Accretion Sweep (15 Nov; journal 113, eighty-seventh pass). Rule fixed 15 Nov before the runs.
rem  Why: STATE section 0, order item 2, the volume. Since 'cover' the continents lose 3.9% of their volume over the
rem  run (band >= -3%). Their losses are Earth's per km of trench; what they lack is Earth's largest addition, island
rem  arcs that reach them (0.00 against ~1.5 km3/yr, Scholl & von Huene 2009). B-d2's parts are built (27-28 Oct,
rem  journals 57-58) and were rejected on 29 Oct only by T2-T4, vetoes then and scored since 4 Nov (journal 59):
rem   trench_polarity 'pair'  one polarity per plate pair, as the overlap resolution consumes crust (the node rule puts
rem                           arcs on the downgoing plate at 22% of oceanic arc nodes): a correctness fix
rem   trench_flux 'consumed'  Earth's trench fluxes are per m2 of ocean consumed; the band law counts node-typed trench
rem                           (1.7-2.7x the consumed area for ocean arcs)
rem   arc_accretion True      consumed arc crust docks on the overriding plate from 15 km (Cloos 1993), else its upper
rem                           3.5 km (Tetreault & Buiter 2014)
rem   arc_H_ocean_m           300 (today's) or 1,000 (Earth's arc rate per m2 consumed: Dimalanta et al. 2002's
rem                           30-95 km3/Myr per km over Bird 2003's intra-oceanic closing rate; journal 57)
rem  Arms (4 x 6 seeds = 24 runs); folders trench_polarity-..._trench_flux-..._arc_accretion-..._arc_H_ocean_m-...:
rem   A node / band / False / 300     today's PARAMS (control; must equal sediment_entrainment-cover of
rem                                   sweep-20261115-141434: identity.txt)
rem   B pair / band / False / 300     the polarity fix alone
rem   C pair / consumed / True / 300  B-d2's three parts at today's arc rate (journal 112's arm): read only
rem   D pair / consumed / True / 1000 the same at Earth's arc rate per m2 consumed
rem  Rule. Vetoes on B and D: V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers close
rem   (crust_budget residuals of the control's size); V4 plate_pieces_review T1 detached area <= 0.01 on every seed;
rem   V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption by mechanism, the first that passes every veto: D > B. D carries Earth's measured flux per m2 consumed
rem   with docking and the polarity fix; B the correctness fix alone. None passes: the control kept, the failing veto is
rem   the next thing to work on. C is read only, not adoptable: under 'consumed' H is a flux per m2 consumed, and 300 m
rem   is no measurement of it (it stood beside the band law's 1.7-2.7x over-count); C shows what the accounting and
rem   docking do at today's H (pre-check: arc growth 3.8 km3/Myr per intra-oceanic km, a quarter of today's).
rem  Full-length cloud pre-check, seed 1017 (unused), control / C / D, read only, not pooled (journal 113 section 3):
rem   bounded on both docking arms (detached area max 0.0020 / 0.0022 / 0.0010; highest node 5,871 / 6,177 / 6,477 m);
rem   D passes A1-A3 (33.2 km3/Myr per km), routed 0.95, arcs reaching the continents over the last 40 Myr 0.00 /
rem   0.18 / 0.62 km3/yr (Earth ~1.5), intra-oceanic trench 25.8k / 30.4k / 17.1k km; C lengthens the trench (129.5k
rem   against 89.2k km at 0 Myr, 34 plates against 21) as docking did on 27 Oct, D does not (92.5k, 17). Watch V2' on
rem   D: its young-ocean gate pooled -317 m on seed 1017 (control -196): the offset scales with the sediment within
rem   ~140 km (flexure under loaded margins read by the gate's local Airy; gate_cause_1017.txt), the watch-point of
rem   STATE section 6. If D fails V2' that way, B can still be adopted and a gate with the model's flexure is next.
rem  Expected (not a band): arcs reaching the continents up, continental area and 20-30 km crust up, > 45 km up
rem   (scored); volume up on the six seeds' mean (seed 1017 alone: -1.6 / -3.2 / -2.8%, trajectories apart).
rem  Scored, reported (no veto): crust_budget (A1-A3, routed share, arcs turned continental + docked, volume),
rem   arc_dock_probe, three_factor, basin, continent, hyps, column_split per arm, thick_tail, lith_crust2,
rem   sediment_where, sediment_by_height; land fraction and p25 / p50 / p90, crust bands, K1-K4, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
set P83=%CAL%\precheck-p83-sediment
set P84=%CAL%\precheck-p84-thicktail
set A=trench_polarity-node_trench_flux-band_arc_accretion-False_arc_H_ocean_m-300
set B=trench_polarity-pair_trench_flux-band_arc_accretion-False_arc_H_ocean_m-300
set C=trench_polarity-pair_trench_flux-consumed_arc_accretion-True_arc_H_ocean_m-300
set D=trench_polarity-pair_trench_flux-consumed_arc_accretion-True_arc_H_ocean_m-1000
%PY% -u -m scripts.run_sweep --name arc-accretion --grid trench_polarity,trench_flux,arc_accretion,arc_H_ocean_m=node:band:False:300,pair:band:False:300,pair:consumed:True:300,pair:consumed:True:1000 --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261115-141434\sediment_entrainment-cover || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% -u -m planet_reviews.arc_dock_probe "%SW%" --control %A% > "%SW%\reviews\arc_dock_probe.txt" 2>&1 || goto fail
%PY% %P80%\lith_crust2.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
for %%a in (A B C D) do call :perarm %%a || goto fail
%PY% %P83%\sediment_where.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" "%SW%\%D%" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:perarm
call set "ARM=%%%1%%"
%PY% %P82%\three_factor.py "%SW%\%ARM%" > "%SW%\reviews\three_factor_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\basin.py "%SW%\%ARM%" > "%SW%\reviews\basin_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\continent.py "%SW%\%ARM%" > "%SW%\reviews\continent_%1.txt" 2>&1 || exit /b 1
%PY% %P82%\hyps.py "%SW%\%ARM%" > "%SW%\reviews\hyps_%1.txt" 2>&1 || exit /b 1
%PY% %P75%\column_split.py "%SW%\%ARM%" > "%SW%\reviews\column_split_%1.txt" 2>&1 || exit /b 1
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
