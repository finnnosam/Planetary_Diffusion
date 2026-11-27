### Borrowed control (21 Nov, after pass 100): the sweeps stop re-running the control; one seed is re-run as the identity check, the other five are taken from the previous adopted arm (user, 21 Nov)

#### 1. Why

Each sweep's control is the previous sweep's adopted arm (the `--identity` reference), and runs on the user's machine
are deterministic, so the control's six runs repeat data already on disk. In the Placed Field Sweep they were 6 of 12
runs and the longest (27–30 min against B's 19–24): with 8 workers they pushed the sweep into a second wave.

The comparison still needs six control seeds: seed-to-seed spread is large (p10 on the datum, sd 238 m between seeds,
more than B's whole effect), and V3 and V5 read the control seed by seed. What it does not need is to *run* them.

The user asked (21 Nov) whether one seed would do; agreed with the reasoning above, and asked for it to be built and
written into the working rules.

#### 2. What was built (no model change)

- `scripts.sweep_plates --borrow NAME=FOLDER [--check-seeds S,...]`: combination NAME runs only its check seeds
  (default the first seed), scheduled first. After the runs each check seed is compared with FOLDER's run of that
  seed by `tools/compare_runs.py` (every shared snapshot field and stats value; `<NAME>/borrowed.txt`).
  - Identical on every check seed: the other seeds are placed in `<NAME>/seed<s>` from FOLDER: snapshots (`.npz`,
    `-plates.json`) and `stats.json` as hard links (files a run writes once and nothing rewrites: no extra disk, and
    deleting either sweep leaves the other whole), run logs and pictures copied, a copy wherever a link fails. They
    are scored like run seeds, so `results.json` and `summary.txt` carry six control seeds.
  - Different on any: the sweep runs the other seeds itself (the old behaviour), so a failed check never stops a
    sweep; the progress line logs it.
  - `<NAME>/borrowed.json` and `sweep.json` (`borrow`) record the source, the check seeds, the verdict, the borrowed
    seeds and the linked / copied file counts. A source missing a seed, or run for other years or another snapshot
    cadence, is refused before anything runs.
- `scripts.run_sweep --borrow-control [--check-seeds ...]` (needs `--control` and `--identity`): passes the borrow to
  `sweep_plates` with the identity reference as the source; the time prediction counts only the runs made, the
  reviews' prediction all runs. Every review and every `.cmd` probe reads the control as an ordinary arm; the
  `identity` review still compares the whole control with the reference (borrowed seeds compare equal by
  construction; the check seed is the real test).
- `tests/test_borrow_control.py` (6, in `PLANET_TESTS`): the job plan, the links and copies, the comparison's
  verdict, the source checks, and `sweep_plates.main` end to end with the runs and Earth's scores stubbed, on a match
  (seeds borrowed, linked) and on a mismatch (the other seeds run).

What a one-seed check can miss: a change that acts only on a rare event (a particular merge or split) that the check
seed does not meet. Each change to the default path also gets the 70-step identity check (n = 48) in the cloud.

Time: a two-arm sweep goes from 12 runs to 7 (one wave on 8 workers), roughly 20–25 minutes shorter; a four-arm sweep
from 24 to 19.

#### 3. The working rule (`docs/next-chat-prompt.md`, Runs and identity)

Every sweep `.cmd` passes `--control <arm> --identity <the previous adopted arm> --borrow-control`. The comparison
stays on six seeds. Disk: hard links share space, so deleting an old reference sweep frees nothing while a later
sweep links its files.

#### 4. Files

Uncommitted on the device: `scripts/sweep_plates.py`, `scripts/run_sweep.py`, `tests/test_borrow_control.py`.
Suggested commit subject: "Borrowed control for sweeps". Docs: this entry, `journal/README.md`, `STATE.md`,
`decisions.md`, `docs/next-chat-prompt.md`, `docs/project-map.md`.
