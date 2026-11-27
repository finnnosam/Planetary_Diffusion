@echo off
rem Lithosphere Frame Sweep (9 Nov; journal 104, eightieth pass). Rule fixed 9 Nov before the runs.
rem  Why: the land's height is the largest failure (land p50 1,165 m against 384; STATE section 0) and the column is
rem  its largest share. C pairs two frames (journal 102 section 4): Artemieva's thermal-age L (93.6 + 0.04 t km) with
rem  Lamb et al. 2020's density excess (30 kg/m3), which Lamb derived for seismic L. lith_frame 'lamb' puts C in one
rem  measured frame: L by tectonothermal age from Hoggard et al. 2020's seismic lithosphere (SL2013sv, data\Hoggard2020)
rem  on Earth's own continent (CRUST1.0 class ages, crust <= 45 km: the planet thickens collided lithosphere itself),
rem  area means per class made monotone (cooling never thins lithosphere): 94 km at 33 Ma, 112 at 100 Ma, 114 to
rem  770 Ma, 165 at 1.5 Ga, 178 for cratons; density excess from Lamb's k (10.753 m per km x 3300: 35.5 kg/m3);
rem  reference 130.68 km and top 177.74 km (thickening bound) from the same table. No fitted number.
rem  Found while building (journal 104): the age law carries less than journal 102 expected. On Earth's ages its
rem  young-old contrast is 65 km (Artemieva's 82); Lamb's 1,020 m between thin-crust and thick-crust regions lives in
rem  L's link to the crust, not to age (Earth's seismic L by crust bin, nodes: 97 / 98 / 107 / 143 / 158 / 157 / 197 km).
rem  Static estimate on the adopted arm's 0 Myr states (static_frame_thicken.txt): land p50 +38 m, p90 +87, the
rem  thick tail +80-125 m (the lower top, 178 against 237.6 km), interiors 35-45 km +28-52 m mean. Built and swept
rem  anyway (build and sweep, don't re-aim); the L-T link is the next arm (STATE section 0).
rem  Cloud pre-check, seed 1011 (unused), 250 Myr, read only, not pooled: both runs bounded (L <= 177.7 km, highest
rem  node 5,778 m); control -> lamb land p50 1,355 -> 1,229 m, sea 97 m higher, 40-45 km crust -235 m, dL/dT 1.66 ->
rem  1.20 km per km (own frames): one seed and a trajectory, against the static +38 m. Rule unchanged.
rem  Expected (written before the runs, not a band): land p50 0 to +80 m, the thick tail higher, C1 a little weaker,
rem  dL/dT in the run's own frame lower than 1.21 (the top at 178 km), L by crust bin flatter; the sweep reads what the
rem  trajectory adds.
rem  Arms (2 x 6 seeds = 12 runs, about 40 min); folders lith_frame-...:
rem   artemieva  today's PARAMS (control; must equal collision_lithosphere-thicken of sweep-20261109-204334: identity.txt)
rem   lamb       the candidate
rem  Rule. Vetoes, every one on its arm: V1 pytest passes (incl. test_lith_frame); V2' elevation gates +-300 m on 6 of
rem   6 seeds; V3 ledgers close; V4 plate_pieces_review T1 detached area <= 0.01 on every seed; V5 on every seed the
rem   highest node over the run <= 8,850 m or no higher than the control's there.
rem  Adoption: lith_frame 'lamb' if its arm passes every veto (one measured frame for L and its density excess
rem   replaces two joined frames: more realistic, whatever a score does; a worse score becomes the next thing to work on).
rem  Scored, reported (no veto): lith_crust2 (each arm in its own frame, beside Lamb's and Hoggard's node references)
rem   and column_split per arm, land fraction and p25 / p50 / p90, sea, C1, crust bands, export, volume, K1-K4, edge,
rem   T2-T5, plates.
cd /d "%~dp0"
set PY=".venv\Scripts\python.exe"
set CAL=outputs\planet-calibration
set P75=%CAL%\precheck-p75-stepback
set P80=%CAL%\precheck-p80-hoggard
%PY% -u -m scripts.run_sweep --name lith-frame --grid lith_frame=artemieva,lamb --seeds 1001,1002,1003,1004,1005,1006 --control lith_frame-artemieva --identity %CAL%\sweep-20261109-204334\collision_lithosphere-thicken || goto fail
for /f "delims=" %%d in ('dir /b /ad /o-n %CAL%\sweep-*') do (set "SW=%CAL%\%%d" & goto found)
:found
echo Probes on %SW%
%PY% %P80%\lith_crust2.py "%SW%\lith_frame-artemieva" "%SW%\lith_frame-lamb" > "%SW%\reviews\lith_crust.txt" 2>&1 || goto fail
%PY% %P75%\column_split.py "%SW%\lith_frame-artemieva" > "%SW%\reviews\column_split_artemieva.txt" 2>&1 || goto fail
%PY% %P75%\column_split.py "%SW%\lith_frame-lamb" > "%SW%\reviews\column_split_lamb.txt" 2>&1 || goto fail
echo Done: %SW%
pause
exit /b 0
:fail
echo FAILED, see the messages above
pause
exit /b 1
