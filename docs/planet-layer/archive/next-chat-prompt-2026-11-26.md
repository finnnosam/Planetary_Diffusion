# Next chat: planet layer

Project: E:\projects\Planetary_Diffusion (git: github.com/finnnosam/Planetary_Diffusion, **public**; Windows, PowerShell,
`.venv` Python 3.12, numba, RTX 3070 8 GB, 16 threads). Split from E:\projects\Infinite_Diffusion on 29 Oct
(journal `60_2026-10-29_repo-split.md`, `61_2026-10-29_packages.md`): that repo keeps the original generator, is not
installed here, and gets no planet work. The migration is complete (30 Oct). Code is one package per stage under `packages/` (`planet_grid`,
`planet_tectonics`, `planet_landscape`, `planet_meso`, `planet_earth`, `planet_reviews`, `planet_climate`); runs in
`scripts/`; stage interfaces in `docs/contracts/` (public: change the contract before a stage's outputs).

**Read first, in this order:**
1. `docs/planet-layer/STATE.md`: **§0 Direction first** (the target, its measured split, the order of candidates, what
   not to do), then what is running, the model as it stands (including the PARAMS-vs-.cmd warning), the scorecard
   against Earth, and the open problems in order. **This says what to do next.** Its order comes from the
   independent review of 22 Nov (`docs/review/2026-11-review.md`, by Astra) and the user's decisions on it:
   `docs/planet-layer/journal/141_2026-11-22_astra-review.md` summarises the review, its three runs and those
   decisions. Read journal 141 next (the review itself only where journal 141 points to it), then the latest passes'
   entries that STATE names (journals 142–148: the age fix, the shoreline trace, the return transfer refuted, the floor
   trace and the remap, the meso on Earth input, the remap's defect and fix, the Remap Sweep and the rift's opened area).
2. `docs/planet-layer/decisions.md`: skim it so you don't re-propose something already rejected (older rows in
   `decisions-archive-oct.md` and `decisions-archive-nov1.md` … `decisions-archive-nov8.md`).
3. `docs/project-map.md` (code, scripts, data, commands) and `docs/planet-parameter-audit.md` (which values are
   measured, assumed or fitted).
4. Journal entries (`docs/planet-layer/journal/`) only when STATE or decisions points to one. Don't read the journal
   front to back. The old single plan doc (`planet-layer-plan.md`) was deleted on 30 Oct: it was exactly
   `design-2026-10-17.md` + journal 01–58 (also in Infinite_Diffusion's tag `pre-split-2026-10-29`). Code comments and
   `.cmd` headers citing "plan doc, Nth pass" mean `journal/*_pass-NN.md`.

## How to keep the docs (every pass)

- **Git:** the GitHub repo is public and for people who download the software (user, 29 Oct). Published: code,
  tests, README, `docs/contracts/`, `docs/planet-ideal-flow.md`, `docs/planet-parameter-audit.md`, the research survey. **Not published**
  (`.gitignore`; on disk and in the Project only): this file, STATE, decisions, journal, project map, plan doc, `.cmd`
  files, `migration/`, `tools/device_*`. So the Project mirror is these notes' only backup: keep it current. Claude does not
  commit (user, 30 Oct: committing is the user's, as it always was). At the end of a pass, leave the code
  changes uncommitted on the device and list them in the journal and the reply: the changed files and a suggested
  commit subject of a few words, such as "Slab initiation rule" (the reasoning lives in the journal). Read-only git in
  the device shell (`git status`, `git diff`, `git log`) is fine; don't ask for delete permission for git. Pushing is
  the user's too.
- **Journal:** write each pass as a new file `docs/planet-layer/journal/NN_YYYY-MM-DD_pass-KK.md` (next: `149_…_pass-111`;
  non-pass entries get a short slug), with a `### ` heading as before. Add its row to `journal/README.md`. Never edit
  an old entry. Correct it in a later one.
- **STATE.md:** rewrite it so it is true after the pass: section 1 (running / to run), the model table, the scorecard
  (state the source sweep), and the open problems. Keep it under ~200 lines.
- **decisions.md:** add one row per adopted / rejected / kept setting, fix or rule, in the same pass. When it passes ~35 KB,
  move its oldest rows, unchanged, to a new `decisions-archive-<period>.md` (each part must fit one read to be mirrored),
  and mirror the new archive as `claude/planet-layer-decisions-archive-<period>.md`.
- **This file:** change it only when a working rule changes.
- **Audit and project map:** update them when values, scripts or `.cmd` files change.
- **Project mirror** (claude.ai Project "Planetary Diffusion"): after a pass, upload `STATE.md` →
  `claude/planet-layer-state-<k>.md` (in parts), `decisions.md` → `claude/planet-layer-decisions.md` (its archives →
  `claude/planet-layer-decisions-archive-<period>.md`, uploaded once), this file →
  `claude/next-chat-prompt.md`, the audit → `claude/planet-parameter-audit-<k>.md` and the project map →
  `claude/project-map.md` if changed, plus the new journal file → `claude/planet-layer-journal/<file name>`.
  Use `project_write` with `local_path` (stage from the device, copy into the cloud workspace). Since 16 Nov the tool
  refuses to replace a doc that this chat has not `project_read` (or written before), and a `local_path` file must
  have been read or written in full in this chat (pass 89): `project_read` the old copy first, and make the local file
  by editing a staged copy you read in full (or write it). Don't read the big Project docs for any other reason.
  **Big docs go up in parts** (user, 19 Nov: "permission to split stuff up"): the Read tool cuts lines over ~2,000
  characters and very large files, so a doc that cannot be read in full in one Read cannot be mirrored. `python3
  tools/device_runs/mirror_split.py <doc>` writes `outputs/mirror/<stem>-<k>.md` (≤ 30 KB, split at `## ` headings; the file
  on disk stays whole) and refuses lines over 1,900 characters: break a long table row into `| ↳ (continued) | | … |`
  rows first. The audit is mirrored this way since 19 Nov (two parts; the old single copy was deleted), the project
  map since 20 Nov (it passed 45 KB; parts `claude/project-map-<k>.md`, the single copy deleted), and STATE since
  26 Nov (pass 110: past one read at 57 KB; parts `claude/planet-layer-state-<k>.md`, the single copy deleted). Upload
  only the parts whose text changed (the first line names the source path, so compare from the second line).
- **New chat:** tell the user when it is a good time to start one, once the docs and the Project copies are current.
  Do this when the chat has grown long or a checkpoint is reached (results read, docs updated, nothing half-built).
  Say it plainly, e.g. "This is a good point to start a new chat; the handoff is up to date."

## Working rules

**Direction and decisions**
- **The direction is fixed in STATE §0** (user, 8 Nov: a change of hands must not drift the work). Follow its order
  and test each candidate against what §0 names (since 22 Nov the geography of STATE §4a: lowland by distance
  inland, four-ring dry area, land–water crossings, one-cell holes; with `basin` and `continent`). Change §0 only by a step-back written
  in the journal and told to the user in the reply (as the ten-pass check); never silently, never by re-aiming a
  built candidate. The user delegated the order within it to Claude (8 Nov); a change of the target itself (the
  land's height) goes to the user.
- **Realism** (user, 23 Oct). Any system may be added, removed or changed if it moves toward realism. Comments may not
  reflect the goal or the code, so decide by the goal, not by a comment. Values come from Earth measurements or from
  mechanism, never from fitting the planet's output to a band. Replace values that only make a band pass.
- **Fitted values only as a last resort** (user, 31 Oct). Compute a new fitted value only when no other good option
  exists (no measurement, no mechanism). Calibrating a mechanism's gain against Earth is a fitted value too.
- **Bring a decision to the user only when it is very important** (user, 26 Oct). That means one of these:
  - it cannot be undone cheaply;
  - it changes the project's goal or direction;
  - it needs something only the user has (data, hardware, money, or a long run on their machine that is not the
    routine sweep).

  Everything else Claude decides by the goal and these rules, and writes the reasoning in the journal: model design,
  which candidate to build next, test bands from Earth data, adopting by a rule fixed before the runs. Running the
  next `.cmd` sweep is routine: just say which one to run.
- Mechanism above ~100 km, learning below. Propose a new model only together with its test set.
- Fix the pass bands and the selection rule before running. Do a total re-evaluation when something looks off.
- **Step back** (user, 31 Oct): when five passes go by without a change to the default path, or every ten passes,
  check the current thread against the terrain's largest failures (STATE's scorecard) before the next pass: is it
  still the most direct route to them, or has it drifted into a sub-question (journal 74: passes 45–52 went three
  steps from the failure they set out to fix)? Write the check in the journal; a change of order goes to the user.
  When stepping back and re-thinking, **run new sweeps, even of changes not expected to help** (user, 4 Nov): a sweep
  gives information no diagnosis can, and at worst it rules something out. Arms still come from a mechanism or a
  measurement (no parameter fishing), the rule is fixed before the runs, and reading an arm is not adopting it.
- **Build and sweep, don't re-aim** (user, 4 Nov; journal 86: passes 53–63 re-aimed eleven times and adopted
  nothing). Every pass ends with something built and swept (a `.cmd` ready to run), adopted, **or an evidence
  checkpoint** (user, 22 Nov, journal 141 §3): a finished causal comparison, read against a question and a reading
  rule written in the journal before the runs, whose result names the next arm. A diagnosis that is not finished
  is not a checkpoint. A candidate re-aimed before it is built is built anyway and the sweep decides; a finding that
  the cause lies elsewhere becomes the next arm, not a reason to drop the current one.
- **"Implemented correctly; geography benefit unresolved"** (user, 22 Nov) is an allowed ending for a candidate: it
  passes the integrity vetoes, its terrain benefit is not shown, so it stays built but off (journal 141 §1: passes
  100, 102 and 104 adopted on integrity while the premise was untested).
- **Vetoes are integrity only; terrain measures are scored** (user, 4 Nov). Vetoes: pytest; the elevation gates
  (±300 m on 6 of 6 seeds); ledgers close; no runaway plates (B-b's T1, detached area ≤ 0.01 per seed); no towers
  (highest node over the run ≤ 8,850 m, or no higher than the control's on that seed). Export, volume, the margins
  band, the hypsometry and every scorecard line are scored and reported: an arm that is more realistic and passes
  the vetoes is adopted even if a score gets worse, and that score becomes the next thing to work on. Adopting a
  change *for its terrain* also needs the confirmation read below.
- Don't judge a law on a transient.
- **Sweeps decide on six seeds (1001–1006), never on one.**
- **Confirmation seeds 2001–2006** (user, 22 Nov, journal 141 §3). 1001–1006 are development seeds: they are read
  again and again, so a design can fit them. Before a change is adopted for its terrain, write in the journal its
  prediction (which measures move, which way, roughly how far), the masks and the time window (endpoint and a
  window declared in advance). Then hand over a `.cmd` that runs the control and the arm on 2001–2006 (borrowing
  does not apply to seeds that were never run). A prediction that fails there is reported and the change is not
  adopted for its terrain. Never use 2001–2006 for design, pre-checks, probes or tuning. Fixes of defects and
  integrity-only changes need no confirmation. Show the endpoint and the per-seed window means together.
- Always compare against `outputs/meso/dataset.npz`, and against ETOPO on the same grid (`planet_earth.reference`).
- Avoid these:
  - learned planners above ~100 km;
  - artificial constraints (locks, frozen sign, projections, mean conservation);
  - scale hacks;
  - over-tuning the test bed;
  - training past the best point;
  - metric gaming.
- Judge terrain on numbers, or on pictures with enough colour levels (`elevation_shading.CMAP` N = 65,536; the old
  256 levels drew false terraces).

**Runs and identity**
- Don't run long generations or trainings yourself. Give simple PowerShell commands with real paths: no
  `<placeholders>`; use `--plates latest` / `--plates sweep`. Sweeps go in a root `.cmd` (the rule in its `rem`
  header) whose body is one `scripts.run_sweep` call: tests, runs and reviews under one progress line with elapsed
  against predicted time (user, 7 Nov: plain progress, not printouts); reviews land in `<sweep>/reviews/`, each run's
  printout in `<arm>/seed<s>/run.log`, step times in `outputs/planet-calibration/timings.jsonl` (the predictions'
  source). A review not in `run_sweep.REVIEWS`, or one that pools the folders it is given (`trench_pairs`, called once
  per arm), goes in the `.cmd` after it, written to a file. Template: `Craton Lithosphere Sweep.cmd`.
- **The control is borrowed, not re-run** (user, 21 Nov; journal 135). Runs are deterministic and each sweep's control
  is the previous sweep's adopted arm, so every sweep `.cmd` passes `--control <arm> --identity <the previous adopted
  arm> --borrow-control` to `run_sweep`: the control runs its first seed only (`--check-seeds` for others), compared
  field by field with the reference (`<control>/borrowed.txt`); identical, its other seeds are hard-linked from the
  reference (no extra disk); different, the sweep runs them itself. The comparison still uses six control seeds, and
  the vetoes and every review read the control as an ordinary arm. Hard links share space: deleting an old reference
  sweep frees nothing while a later sweep links its files. The borrow needs the reference run with numerically
  identical code at the sweeps' `threads` 2 (pass 102: a speed-up checked at `threads` 1 broke the borrow at the 14th
  digit; journal 137 §1). Until pass 104 a run at `threads` 2 did not repeat itself (ducc0's threaded adjoint, the slab
  parcels' coefficients from step 7; journal 140 §1), so no borrow had ever succeeded; the adjoint now runs on one
  thread. Check a change claimed "identical" at `threads` 1 and at `threads` 2, and a run against itself at `threads` 2,
  with `precheck-p104-repeat/rep.py` in the device shell (n = 48, 70 steps, two runs side by side, ~2.5 min); the
  sweep's borrow check seed tests it again on Windows.
- Runs on the user's machine are deterministic (at `threads` 2 since pass 104's fix; `Repeatable Control.cmd` tests it). The cloud workspace and the device shell follow their own
  trajectories (seed 1001 in the cloud: 26–28 plates at 0 Myr vs 20 on Windows). Use them for pre-checks only.
- **Cloud pre-checks run on seeds outside the sweep's six** (user, 31 Oct): 1007 first, then 1008, ... (1033 is the next unused; 1032 was the pass-110 full-length pre-check, 1031 the pass-109 one, 1030 pass 108's, 1029 pass 107's, 1028 the pass-105 trace and its warm state, re-used in passes 106–108; 2001–2006 are reserved for confirmation) (never used
  before 31 Oct). A cloud run of seed 1001 shares its start with Windows' 1001, so it is a correlated sample; an
  unused seed is an independent one, and the pre-checks stop tuning toward one seed. Run its control beside it (the
  old cloud-1001 numbers do not carry over). Report it beside the six-seed means as a seventh read; it is not pooled
  into a rule fixed on six seeds. Paired probes may use any seed (they compare arms on the same states).
- **A new rule that accumulates state** (a field multiplied or added to every step) **gets a full-length cloud
  pre-check** (seed 1007+, 250 Myr) before its sweep is handed over: short tests miss ratchets (journal 102). Read
  its run logs for `RuntimeWarning` beside the vetoes: a warning the control does not have is investigated before the
  sweep is handed over (pass 109: two overflow warnings per run were the only sign of a return that went to inf on one
  Windows seed).
- Runs drift apart from the first steps. Measure an effect smaller than that drift paired, on the same states.
- Step-by-step identity checks need `params={'threads': 1}` (ducc0 threads differ at 1e-12). When editing the default
  path, check that it stays identical step by step (n = 48, ~70 steps). Call `sim.state()` once before the first step
  and after every step, as `simulate_plates` does: `state()` sets `z_last`, which the measured dynamic topography
  reads at the next step (journal 88), so harnesses that call it at another cadence follow another trajectory.
- Sweep folders are named only by the parameters that differ (`MAX_NAME` 80, else `combo-NN`). Shared parameters go in
  `sweep.json`, per-combination ones in `<combo>/params.json`. Runs save `params_effective`. Reviews read older runs
  with `crust_review.LEGACY_DEFAULTS`.
- `planet_meso_check --landscape latest` prefers the newest theta calibration, so pass explicit landscape folders.
- Continental crust in the reviews is label-free: total (crust + sediment) ≥ 20 km on both sides (`crust_review` at
  46 km; `margin_width` on CRUST1.0's 1° cells via `to_crust1_cells`). The flag reading is printed beside it.
- Tests: `pytest -q tests` (all, meso included; Windows). Planet-layer subset: `pytest -q tests\test_shelf.py
  tests\test_plates.py tests\test_sphere_graph.py tests\test_graph_landscape.py tests\test_coupling.py
  tests\test_boundary_lengths.py tests\test_erosion_law_check.py tests\test_crust_reviews.py tests\test_trench_pairs.py
  tests\test_mantle_flow.py tests\test_flow_reference.py tests\test_save_cadence.py tests\test_lithosphere.py
  tests\test_dyn_measured.py tests\test_margin_wedge.py tests\test_collision_spread.py` (the sweeps run the full
  list, `scripts/run_sweep.py` PLANET_TESTS, incl. `test_pelagic`, `test_deep_crust`, `test_collision_lithosphere`, `test_lt_link`, `test_ocean_depth`).
- Sweeps save every 5 Myr over the last 40 (`sweep_plates --save-every`, default 5; pass 48). The landscape stage needs
  every Myr: rerun the chosen seed with `simulate_plates` (default 1).
- `PARAMS` are the working baseline (since 29 Oct): a bare run is the current model. New code reads `data/` and
  `outputs/` through `planet_grid.paths` (`DATA`, `OUTPUTS`), never from paths relative to a file.
- The user's `.venv` has no scipy or scikit-image, so planet-layer code stays numpy / numba / ducc0 (+ matplotlib).

**The user's machine through the device shell** (device_bash; repo at `$HOME/mnt/Planetary_Diffusion`)
- Run scripts from the repo root: `planet_grid.paths` defaults (`data`, `outputs`) are relative to the working
  directory (pass 52).
- Reviews and probes run there with `NUMBA_CACHE_DIR=$HOME/numba-cache PYTHONDONTWRITEBYTECODE=1` and
  `PYTHONPATH=$PWD:$(ls -d $PWD/packages/*/src | tr '\n' ':')` (the packages are not installed there; pytest finds
  them through `conftest.py`). Each call is ≤ 180 s, and background jobs are stopped. It has 2 cores and
  3.9 GB. Save review printouts beside the sweep (`<sweep>/reviews/*.txt`).
- Its own disk (`/sessions`, ~10 GB, the app's `sessiondata.vhdx`) only grows: it was full on 6 Nov (136 MB free)
  and the user reset it (renamed `%LOCALAPPDATA%\Claude-Data\claude_pzs8sxrjxfjjc\vm_bundles\claudevm.bundle`; the app
  rebuilt it). If `df -h $HOME` is near full again, tell the user; meanwhile run probes in the cloud (stage per-arm
  tarballs of ≤ ~40 MB; one 218 MB tar timed out). Its Python packages can vanish between sessions. Reinstall with `pip install --user --break-system-packages numpy
  numba ducc0 pytest h5py`. It has no torch, tkinter, rasterio or pyfastnoiselite, so the meso, `arrangement_metrics`
  and the Earth raster readers don't run there. Without numba the reviews still run on the repo's pure-Python fallback.
  `python3 tools/device_stubs/run.py <test names>` smoke-tests plate code with stubs. That is plumbing only, never
  for judging.
- Full n = 217 runs there use `tools/device_runs/chunk_run.py` + `simpickle.py`. They are pickled between 180 s
  calls, two in parallel take ~12 min per run, and they are for pre-checks only.
- `test_coupling` must run from a copy (`$HOME/work/tc/tests`): the shell cannot delete its old `outputs/test-tmp`.
  The shell cannot delete files at all without the user's permission. For the same reason pytest there needs
  `--basetemp=$HOME/work/<new name>` (`mkdir -p $HOME/work` first) and `-p no:cacheprovider` (pass 74).
- Code changes: make them in a cloud clone (commits there are scratch; right after cloning run `git remote remove origin`:
  the cloud workspace's stop hook asks to commit and push whenever the shell sits in a repo with a remote, and
  pushing is the user's, user 7 Nov) and copy the changed files to the device
  with `device_commit_files`; check each file's md5 on both sides. A commit onto an existing device file can report
  "written" and leave the old content (31 Oct): commit to a new path (`outputs/<name>`) and `mv` it into place in
  the device shell. Never run git commands that write in the device's
  repo. Small edits to ignored files (docs) in place:
- Edit files in place with short Python read-modify-write scripts that assert each old string matches once. Check the
  md5 before and after.
- Long shell commands fail (`spawn E2BIG`). Write big files in the cloud workspace and copy them with
  `device_commit_files` from `/mnt/user-data/outputs/<fresh folder>/<new name>`. Use a new file name for every sync
  (a re-commit from the same staged path can write old content). Check the device md5 before overwriting and verify
  it after.
- `pkill -f` / `pgrep -f` patterns can match your own command: pkill kills your shell, and a wait loop on pgrep
  never ends (both happened again in pass 65, in the cloud shell too). Kill or wait by PID from
  `ps aux | grep "simulate_[p]lates"`.

**Cloud workspace**
- It can do real runs: `pip install --break-system-packages numpy numba ducc0 pytest`, then the code (a clone of
  the GitHub repo, or `packages/`, `scripts/`, `tests/`, `tools/`, `conftest.py` staged as tarballs) with
  data/crust1, data/tectonicplates-master and data/mantle (the mantle flow reads `visc.sh08`; 1.6 MB). Set `PYTHONPATH` to the repo and every `packages/*/src`. There is no ETOPO there, so `elevation_check`'s command does
  not run; its gates do (`elevation_check.snapshot_tests`, pooled as `check()`: `precheck-p94-level/precheck_gate.py`).
  Staging `outputs/planet-reference/` `earth-n217.npz` / `.json`, `earth-plate-speeds-nnr.json`, `earth-plates-n217.npz`,
  `earth-seafloor-age-n217.npz`, `earth-trenches.json`, `globsed-n217.npz` (~2.5 MB) lets `crust_budget` (V3) and
  `plate_pieces_review` (V4) run there too (pass 94). Two 250 Myr n = 217 runs in
  parallel with `--param threads=1` take ~25–60 min (Bash `run_in_background` with `timeout` 7200000: the default
  limit, 30 min, stops the runs; pass 87 lost a first try at −15 Myr that way).
- Neither shell reaches most data servers. The user downloads data in a browser.
- pygplates installs from PyPI in the cloud (`pip install --break-system-packages pygplates`; 1.0.0). Earth's plate
  model and paleo-age grids are in `data/Muller2019/`; age grids are sampled on the device (h5py; E3, pass 49).
- BurnMan (SLB2011 mineral physics) installs in the cloud with `pip install --break-system-packages --no-deps burnman`
  plus scipy (its numba pin fails to build; it warns about numpy 2 and runs). Cloud only; the package carries tables.
- Compiling external code (e.g. HC) is refused by the permission classifier. Don't work around it: write the
  reference from scratch, sharing no code or convention with what it checks (pass 48).
