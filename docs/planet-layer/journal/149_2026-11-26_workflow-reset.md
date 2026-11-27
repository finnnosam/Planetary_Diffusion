### Workflow reset (26 Nov, after pass 110): Claude Code on the Windows machine, the notes in git, a one-page handoff, and the target restated as what the map shows

No model change and no runs. The user asked for a review of whether the instructions and workflow, including anything
the user had asked for, were getting in the way of the work, keeping only the project's goal fixed (a realistic
Earth-like planet generator; the plate-tectonics stage first). This entry records what was found, what the user
decided, and what was changed.

#### 1. What was getting in the way

- **The agent could not run its own experiments.** Passes 1–110 ran from a cloud chat linked to the machine through a
  small Linux VM (2 cores, 3.9 GB, 180 s per call, no deletes, no git writes). The VM and the cloud workspace follow
  other trajectories than Windows, so they served only for pre-checks; every sweep went to the user as a `.cmd`, and a
  new chat read it. About 55 of the 246 lines of `next-chat-prompt.md` were workarounds for that link (copies to fresh
  paths with md5 checks, `E2BIG`, the commit bug, pytest's `--basetemp`, pkill matching itself). The user had
  understood that Claude Code could not run on Windows; it can (the native installer).
- **Bookkeeping cost about as much as the work.** The reading list before the first action was ~340 KB (~90k tokens).
  STATE had grown to 482 lines and 57 KB against its own ~200-line rule. The journal is 1.6 MB in 149 entries; its index
  62 KB; the decision log 475 rows in ten files, all to be skimmed. Each pass updated six or seven documents and mirrored
  them to the claude.ai Project in under-30 KB parts. The root cause: the notes were kept out of the public repo, which
  made the Project their only backup.
- **The target's measurement code was not versioned.** `geography.py`, the measure §0 judged by since pass 107, lived
  in `outputs/planet-calibration/precheck-p107-transfer/`, as did the probes the sweeps call from thirteen
  `precheck-pNN` folders. `outputs/` is ignored by git and the Project mirror held documents only: no backup, no tests.
- **The rules had built up one incident at a time** (17 dated user rules, a 13-item "not to do", codes such as V1–V5,
  B-d2, T2–T5, K1, and "C" for different arms in different passes). Each is sound alone; together a fresh agent spends
  much of its effort complying and decoding.
- **Self-run step-backs confirmed the current thread.** Journals 117 and 123 found the inland-lowland problem, but the
  passes that followed did not follow it (the review's finding); journals 147 §5 and 148 §6 found "the order stays".
  The step-back that changed course was the independent review (journal 141), which read code and data first.
- **Confirmation cost a round trip that was not needed.** Journal 145 §6 froze its prediction before the Remap Sweep, so
  2001–2006 could have run in the same sweep.

#### 2. The user's decisions (26 Nov)

1. Work moves to Claude Code on the Windows machine.
2. The public repo may hold more than what a user needs to run the code: the notes, the sweep `.cmd` files and the
   probe scripts go in; runs and arrays stay out (no 100 GB of experiments).
3. The target. Asked, the user was not sure what the continents' level on the datum signifies. It is the height of
   continental crust against the fixed 0 m, not against the model's own sea: a diagnostic that separates continents
   standing too low from a sea standing too high, which nobody sees on the map. Claude recommended making the target
   what the map shows (land fraction; where the land is; how high it is, relative to the sea), with the datum level as a
   diagnostic, the target deciding when they disagree. The user answered "set up everything", taken as yes; it is the
   user's to revert. This replaces the target chosen on 19 Nov (journal 122 §6).
4. Six seeds stay: the user considers six enough; to be settled more rigorously some time.
5. No regular end-to-end renders now: the plate stage needs more work first, and a render would only amplify its
   problems.
6. The ban on fitted values stays, for the user's reason: fitted values interact, so every new one throws off those
   fitted before it.

#### 3. What was changed

- `CLAUDE.md` at the root (read automatically by Claude Code) replaces `docs/next-chat-prompt.md` (moved to
  `archive/next-chat-prompt-2026-11-26.md`). One page: the goal, the target, how work is judged, what goes to the user,
  how to run things on Windows, how to keep the notes. New in it: Claude runs routine sweeps itself and says when it
  starts one; confirmation seeds run in the same sweep when the prediction is written first; an independent review by a
  fresh subagent every ~15 passes replaces the self step-back; Claude commits code and notes at the end of a pass (the
  user had kept commits; pushing stays the user's); STATE's limit is ~120 lines.
- `STATE.md` rewritten to one page; the old one moved unchanged to `archive/STATE-2026-11-26-pass110.md`.
- `gotchas.md`: the traps, collected unchanged from the old STATE §6 (watch-points), its queued corrections, and the
  working rules on runs, identity and borrowing that still apply.
- `.gitignore`: the notes and the sweep `.cmd` files are tracked (the two migration `.cmd`s and `migration/` stay out);
  under `outputs/` only the `.py` and `.md` files of `planet-calibration/precheck-*/` and the `.py`, `.md` and `.ps1`
  files of `review-astra/` (480 new files, 2.8 MB). `tools/device_runs/` and `tools/device_stubs/` stay out.
- `planet_reviews.geography` (new, from `precheck-p107-transfer/geography.py`, same code; it also takes a sweep folder
  with `--control`): read on the Remap Sweep it matches that sweep's `reviews/geography.json` to 3e-14.
  `tests/test_geography.py` (2 tests: a known map, a sweep folder), added to `run_sweep.PLANET_TESTS`; the review added
  to `run_sweep.REVIEWS`. The old script is a shim calling the module, so the older `.cmd` files still run.
- `docs/project-map.md`: the working rules' new home, the geography review, and the cloud-session section replaced.
- The Project mirror stops. The Project's copies (`claude/planet-layer-*`, `claude/project-map-*`,
  `claude/planet-parameter-audit-*`, `claude/next-chat-prompt.md`) can be deleted once this is committed and pushed.

Suggested commit: "Workflow reset: notes in git, geography review" (changed: `.gitignore`, `CLAUDE.md`,
`scripts/run_sweep.py`, `packages/planet-reviews/src/planet_reviews/geography.py`, `tests/test_geography.py`, and the
newly tracked notes, `.cmd` files and probes).
