# Planetary Diffusion: working notes for Claude

The goal is a generator of realistic Earth-like planets. The current focus is the plate-tectonics stage
(`packages/planet-tectonics`, 46 km grid, 250 Myr). Later stages (landscape, meso, climate) wait until it is good.
Claude works here as Claude Code on the owner's Windows machine (16 threads, 32 GB, RTX 3070), in this repo.

## Start of a session

Every session starts fresh. What carries over is only what is in files: this `CLAUDE.md` and the notes it points to,
all in git. Claude Code's automatic memory is off for this project (`.claude/settings.json`), so nothing is kept
anywhere else; if something should last, it goes into the notes. If the user finds a session acting on something
these files don't say, check `/context all` for other memory files (a `CLAUDE.md` in a parent folder or in
`~/.claude/`).

1. Read `docs/planet-layer/STATE.md` (one page: target, scorecard, what is running, what is next).
2. Read only the journal entries and files STATE points to. Search `docs/planet-layer/decisions*.md` with grep
   before proposing a setting, so you don't re-propose something already rejected.
3. `docs/planet-layer/gotchas.md` collects the traps found so far (numerics, borrowing, reviews, reading Earth data).
   Read the parts that touch what you are about to change.
4. Then do the next pass. The user's opening message is usually just "Continue" (or "next pass", "go on"): that
   means do the next pass, numbered one after the last journal entry, on what STATE's "Right now" names as next. If
   a sweep is still running, check on it first. Don't ask the user what to work on; STATE says.

Other references: `docs/project-map.md` (code, scripts, data), `docs/planet-parameter-audit.md` (where each value
comes from), `docs/planet-ideal-flow.md` (the design of every stage), `docs/contracts/` (what each stage hands on).

## The target

What a person would see on the map, against Earth on the same grid (ETOPO, ice regions left out), relative to sea level:

- how much land there is (land fraction);
- where it is: broad inland lowlands and dry interiors, coastline roughness (land-water crossings, one-cell holes and
  islands);
- how high it is: the land's height distribution, mountains and high plateaus.

These are measured by `planet_reviews.geography` and the hypsometry reviews, which every sweep runs. The continents'
level on the datum (height against the fixed 0 m, not against the model's sea) is a **diagnostic**: it tells whether
missing land comes from continents standing too low or from the sea standing too high. When a change improves the
target but worsens the datum level, the target decides; the regression becomes the next problem.
(User, 26 Nov, journal 149. Changing the target itself is the user's decision.)

## How the work is judged

- **Realism from mechanism or measurement.** Any part of the model may be added, changed or removed if it moves toward
  realism. Values come from Earth measurements or from a mechanism. Fitted values are a last resort, and so is
  calibrating a mechanism's gain against Earth: fitted values interact, so each new one throws off the ones fitted
  before it. Never fit a value so that an output lands in a band.
- **Fix the reading rule before the runs.** Write the question, the pass bands and the rule that will decide in the
  journal before starting a sweep or a probe. Reading an arm is not adopting it.
- **Vetoes are integrity only:** pytest; the elevation gates (±300 m on 6 of 6 seeds); the ledgers close (V3); no
  runaway plates (detached area ≤ 0.01 per seed); no towers (highest node ≤ 8,850 m, or no higher than the control's).
  Everything else is scored and reported. An arm that is more realistic and passes the vetoes is adopted even if a
  score gets worse; that score becomes the next thing to work on.
- **Seeds.** Sweeps decide on development seeds 1001–1006, never on one seed. Effects smaller than the drift between
  runs are measured paired, on the same states. Don't pool one run's epochs as independent samples.
- **Confirmation seeds 2001–2006** are never used for design, probes or tuning. Before adopting a change for its
  terrain, write its prediction (which measures move, which way, roughly how far; masks; endpoint and window) in the
  journal, then run control and arm on 2001–2006. When the prediction is written before the development sweep, run
  both seed sets in the same sweep and read the development seeds first. A failed prediction is reported and the change
  is not adopted for its terrain. Defect fixes and integrity-only changes need no confirmation.
- **A pass ends** with something built and swept, adopted, or a finished evidence checkpoint (a causal comparison read
  against a question fixed beforehand, whose result names the next arm). "Implemented correctly; geography benefit
  unresolved" is an allowed ending: the option stays built and off. Don't re-aim a built candidate: sweep it.
- **Independent review** every ~15 passes, or when several passes go by without progress on the target: start a fresh
  subagent that reads code and data **before** STATE and the journal (as `docs/review/2026-11-review.md` did), and
  bring its findings to the user. This replaces the old self-run step-back.
- **Avoid:** learned planners above ~100 km (mechanism above ~100 km, learning below); artificial constraints (locks,
  frozen signs, projections, mean conservation); scale hacks; reading a band selected on elevation as the cause of its
  own excess; judging on land median alone or on one seed; a conditional texture score without its qualifying area.
- Judge terrain on numbers, or on pictures with enough colour levels (`elevation_shading.CMAP`).

## What goes to the user

Only decisions that cannot be undone cheaply, change the goal or the target, or need something only the user has (data
downloads, money, a long non-routine run such as 500 Myr). Everything else Claude decides by the goal and these rules
and writes the reasoning in the journal. Routine sweeps, tests and probes Claude runs itself.

## Running things

- Python is `.venv\Scripts\python.exe` (Python 3.12; numpy, numba, ducc0, matplotlib; no scipy or scikit-image, so
  planet code stays on those). Run everything from the repo root: `planet_grid.paths` resolves `data/` and `outputs/`
  from the working directory.
- Tests: `.venv\Scripts\python.exe -m pytest -q tests` (all); sweeps run the planet subset `scripts/run_sweep.py`
  PLANET_TESTS (~6 min).
- **Sweeps:** one root `.cmd` per sweep, its rule in the `rem` header, its body one `scripts.run_sweep` call followed by
  any extra probes (template: `Rift Area Sweep.cmd`). Every sweep passes
  `--control <arm> --identity <the previous adopted arm's folder> --borrow-control`: runs are deterministic, so the
  control re-runs its first seed only, checks it field by field, and links the other seeds. Run a sweep yourself **in a
  tab of the user's Terminal panel** (`run_in_terminal`, cwd the repo root, titled with the sweep's name), so the user
  sees its live progress bar and time left; link the tab when you start it, read it with `read_terminal` to check on
  it, and stop any watcher you start when the sweep ends. The `.cmd` files end in `pause`, so feed them NUL:
  `cmd.exe /c 'call ".\Rift Area Sweep.cmd" < NUL'`. The `.\` is needed where the shell sets
  `NoDefaultCurrentDirectoryInExePath=1` (Claude Code's own shell does; pass 111, journal 150: without it cmd answers
  "is not recognized"). A borrowed sweep takes ~30 min, a full one about an hour, with 8 workers and ~23 GB of memory:
  tell the user when you start one, since the machine is theirs too.
- Results: `outputs/planet-calibration/<sweep>/`, each run's printout in `<arm>/seed<s>/run.log`, reviews in
  `<sweep>/reviews/`. Read `run.log` for `RuntimeWarning` and `nan`.
- **Determinism:** runs repeat exactly on this machine at `threads` 2. A change claimed to leave the default path
  identical is checked step by step at `threads` 1 and 2 (`precheck-p104-repeat/rep.py`, n = 48, 70 steps).
- **Pre-checks:** a new rule that accumulates state gets one full-length run (250 Myr) on an unused seed before its
  sweep (next unused: 1033; record it in STATE). Run its control beside it. It is a seventh read, never pooled.
- Earth references on the same grid are in `outputs/planet-reference/` (`planet_earth.reference`). Data downloads are
  the user's (no data servers are reachable from scripts).

## Keeping the notes (every pass)

- **Journal:** a new file `docs/planet-layer/journal/NN_YYYY-MM-DD_pass-KK.md` (non-pass entries get a short slug),
  with a one-line row in `journal/README.md`. Write in plain sentences; name measures in words, not codes. Never edit
  an old entry; correct it in a later one.
- **STATE.md:** rewrite it so it is true after the pass. Hard limit ~120 lines: target, current model in a few lines,
  scorecard (about a dozen lines, with the sweep it comes from), what is running, open problems in order. History
  belongs in the journal.
- **decisions.md:** one row per adopted / rejected / kept setting or fix. Past ~35 KB, move the oldest rows unchanged to
  a new `decisions-archive-<period>.md`.
- Audit and project map: update them when values, scripts or `.cmd` files change. A probe that becomes a standing
  measure moves into `planet_reviews` with a test and a `run_sweep.REVIEWS` entry.
- **Git:** code, tests, docs, the notes, the sweep `.cmd` files and the probe scripts in
  `outputs/planet-calibration/precheck-*/` are tracked; data, runs, arrays and printouts are not (`.gitignore`). At the
  end of a pass Claude commits code and notes together in one commit (a subject of a few words, such as "Rift area by
  the plates' opened area"; the reasoning lives in the journal). Then push it to `origin main` (user, 26 Nov: every commit is pushed). The
  repo is public, so before pushing check that the commit holds only tracked kinds of file (no data, runs, arrays or
  secrets); if a push fails, say so to the user rather than forcing it.
- Tell the user when it is a good time to start a new session: the notes are current and nothing is half-built.
