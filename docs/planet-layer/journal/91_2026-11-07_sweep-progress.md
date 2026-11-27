### Sweep progress (7 Nov, after the sixty-eighth pass): one progress line instead of printouts; every step's time recorded

The user asked for sweeps to report progress plainly: elapsed against predicted time, a percentage for the whole and
for its steps, simple lines such as "completed run 1, settings 1, seed 1", a progress bar that updates in place in the
console, and every step's time recorded so later predictions stay accurate.

**Built** (output only; the model and its runs are unchanged):
- `scripts/progress.py`: `Tracker` draws one line, redrawn in place with `\r` (plain lines every 5% when the output is
  not a console), e.g.
  `[#########-----------]  45% | 0:27:10 of ~1:00:20 | 2/3 sweep: 9 of 24 runs done, 8 running [41%]`, with event
  lines printed above it. `record()` appends each finished step to `outputs/planet-calibration/timings.jsonl`;
  `predict()` / `predict_run()` take the median of recent matching records (runs by workers and length; reviews per
  run reviewed), else defaults measured on `sweep-20261107-094624`. `schedule_left()` turns per-run times into the
  pool's remaining time (running runs at their own pace, the queue onto the worker free first).
- `scripts/run_sweep.py`: the whole sweep for a `.cmd`: 1 tests (the planet-layer subset `PLANET_TESTS`, its share read
  from pytest's own `[ NN%]`), 2 the runs (`sweep_plates`), 3 the reviews (`REVIEWS`, as the `.cmd` files ran them since
  pass 67, each into `<sweep>/reviews/`). At the end `<sweep>/reviews/timings.txt` lists each step's predicted and
  actual time. Exit status 1 on any failure (failed tests are printed).
- `scripts/sweep_plates.py`: workers report each simulated Myr through a queue and write their printout to
  `<arm>/seed<s>/run.log` (no longer interleaved on the console); a line per finished run ("run 3 of 12 done:
  settings 1 of 2 (lithosphere-none), seed 1003 (3 of 6), 18.2 min"); the score tables go to `<sweep>/summary.txt`;
  `main()` returns the sweep folder. Without a tracker it draws its own line.
- `scripts/simulate_plates.py`: a module hook `PROGRESS` called with the run's time after each recorded step.
- `tests/test_progress.py` (6 tests; in `PLANET_TESTS`): formatting, the pool schedule, records and predictions, both
  drawing modes, the review command lines, and the hook output only (n = 48, 3 Myr, threads 1: stats and every
  snapshot field equal with and without it). With `test_save_cadence` they pass in the cloud and in the device shell.
- `timings.jsonl` on the device seeded from the last three sweeps (84 runs, median ~19 min; the reviews of
  `sweep-20261107-094624` from their file times; pytest 117 s), marked `bootstrap`. Prediction for the Craton
  Lithosphere Sweep: ~42 min (tests ~2, runs ~38 in two waves, reviews ~2).
- `Craton Lithosphere Sweep.cmd`: same rule and arms; its body is one `scripts.run_sweep` call. Cloud smoke run (2
  settings x 2 seeds x 3 Myr, the scorer stubbed: the cloud lacks torch and rasterio): the second run's prediction
  from the first run's records was 37 s against 38 s actual.

The `.cmd` no longer types every review to the console; the reviews are read from `<sweep>/reviews/` as before.

**Code changed** (uncommitted): new `scripts/progress.py`, `scripts/run_sweep.py`, `tests/test_progress.py`; changed
`scripts/sweep_plates.py`, `scripts/simulate_plates.py`. Suggested commit: **"Sweep progress line and timing record"**.
Not published: the `.cmd`, `timings.jsonl`, docs.
