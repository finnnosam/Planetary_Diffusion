@echo off
rem Start Sediment Sweep (16 Nov; journal 115, eighty-ninth pass). Rule fixed 16 Nov before the runs.
rem  Why: STATE section 0, order item 2, the margins and the low continents. The land stands ~0.12 km too high and
rem  the continents' spread on the datum is 109% of that excess, on the low side (a third of the continent >= 30 km
rem  under the datum against Earth's 11%); the low slice lacks -310 m of column in sediment, and land above the sea
rem  holds 88 against Earth's 176 Mkm3. The start carries no sediment at all, while Earth's continents hold 320 Mkm3
rem  and its floor ~190 (journal 108 section 3), much of it older than the run.
rem  start_sediment (built this pass; tests\test_start_sediment.py):
rem   'none'       as before (control)
rem   'continent'  the supercontinent carries the cover of the Earth domain its crust is already mapped onto
rem                (CRUST1.0 sedthk on Precambrian and platform crust > 500 km from Bird's boundaries: p25 / 50 / 75 /
rem                90 = 100 / 600 / 1,810 / 3,800 m bulk, mean 1,331 m; 267 Mkm3 bulk on seed 1018's start), placed by a
rem                Gaussian copula with Earth's rank correlation to the crystalline crust (-0.31) and Earth's spatial
rem                correlation of the cover (0.84 / 0.50 / 0.21 at < 100 / 200-400 / 600-1,000 km); own random stream,
rem                so the start's crust, plates and poles are unchanged
rem   'all'        the same, and the start's floor carries the pelagic rain of its own age (pelagic_m_myr x age:
rem                36 Mkm3 bulk on seed 1018)
rem  Arms (3 x 6 seeds = 18 runs): start_sediment-none (control; must equal arm D of sweep-20261115-213654:
rem   identity.txt), start_sediment-continent, start_sediment-all.
rem  Rule. Vetoes on 'continent' and 'all': V1 pytest passes; V2' elevation gates +-300 m on 6 of 6 seeds; V3 ledgers
rem   close (crust_budget residuals of the control's size); V4 plate_pieces_review T1 detached area <= 0.01 on every
rem   seed; V5 on every seed the highest node over the run <= 8,850 m or no higher than the control's.
rem  Adoption by mechanism, the first that passes every veto: 'all' > 'continent' ('all' gives both domains of the
rem   start what they carry: Earth's measured cover and the model's own floor law). None passes: the control kept, the
rem   failing veto is the next thing to work on.
rem  Full-length cloud pre-check, seed 1018 (unused), control / 'all', read only, not pooled (journal 115 section 4;
rem   precheck-p89-startsed\precheck_read_1018.txt, gate_1018.txt, t3_origin_1018.txt): every veto reading passes on
rem   'all' (gates pooled young-ocean -198 / thin-crust -218 m against the control's -189 / -135; detached area max
rem   0.0014; highest node 6,082 m against 6,305). The start's cover leaves the high land within ~40 Myr (continental
rem   sediment above the sea 234 -> 130 Mkm3 by -210 Myr) and settles on the low continents and the margins. At 0 Myr
rem   (one seed; trajectories apart, 29 against 40 plates): land p50 418 against 563 m (Earth 380), continents >= 30
rem   km on the datum p25 -58 against -188 m (Earth 109), share under the datum 0.31 against 0.38 (0.11), the sea -96
rem   against -128 m; sediment on continents above the sea 106 against 83 Mkm3 (176). Watch V2' thin-crust on 'all'.
rem  Expected (not a band): sediment on the continents up, most of all on the low continents; the continents' volume
rem   (crust + sediment) down against the control as the start's cover is exported (read the crust alone beside it);
rem   the sea up early (floor and drowned cover displace water).
rem  Scored, reported (no veto): three_factor, basin, continent, hyps, column_split per arm, sediment_where,
rem   sediment_by_height, thick_tail, lith_crust2, crust_budget, T3 by the origin of its continental side (t3_origin,
rem   journal 114 section 4); land fraction and p25 / p50 / p90, crust bands, K1-K4, T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
set P82=%CAL%\precheck-p82-stepback
set P83=%CAL%\precheck-p83-sediment
set P84=%CAL%\precheck-p84-thicktail
set P89=%CAL%\precheck-p89-startsed
set A=start_sediment-none
set B=start_sediment-continent
set C=start_sediment-all
%PY% -u -m scripts.run_sweep --name start-sediment --grid start_sediment=none,continent,all --seeds 1001,1002,1003,1004,1005,1006 --control %A% --identity %CAL%\sweep-20261115-213654\trench_polarity-pair_trench_flux-consumed_arc_accretion-True_arc_H_ocean_m-1000 || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P80%\lith_crust2.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P84%\thick_tail.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" > "%SW%\reviews\thick_tail.txt" 2>&1 || goto fail
for %%a in (A B C) do call :perarm %%a || goto fail
%PY% %P83%\sediment_where.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" > "%SW%\reviews\sediment_where.txt" 2>&1 || goto fail
%PY% %P83%\sediment_by_height.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" > "%SW%\reviews\sediment_by_height.txt" 2>&1 || goto fail
%PY% %P89%\t3_origin.py "%SW%\%A%" "%SW%\%B%" "%SW%\%C%" > "%SW%\reviews\t3_origin.txt" 2>&1 || goto fail
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
