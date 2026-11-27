### Migration complete (30 Oct): the final Windows check passes, seed 1001 IDENTICAL on the plain defaults; Infinite_Diffusion is the generator alone and its tests pass

Follows [61](61_2026-10-29_packages.md). No model change.

**Phase B's final check** (`Finish Migration.cmd`, Windows, 29 Oct evening):
- `pip install -r requirements-dev.txt`; `spherical-terrain-diffusion` uninstalled from this `.venv`.
- `migration/import_check.py`: 77 package modules import, 486 import statements in scripts, tools and tests resolve,
  0 problems, and `planet_diffusion` is not importable.
- All tests: 146 passed (meso included).
- `python -m planet_earth.manifest --record`: all required files present; `data_checksums.json` committed (`c0371ac`).
- Identity: `sweep_plates --rule shelf --seeds 1001` with no `--grid` (`sweep-20261029-215013/defaults`) against
  `sweep-20261028-200644/trench_polarity-node_trench_flux-band_arc_H_ocean_m-300_arc_accretion-False/seed1001`, with
  `--ignore extension_age_myr,thickening_age_myr`: IDENTICAL, 62 snapshots × 16 fields. `stats.json` differs only by
  removal: the extension statistics (`extension_starts`, `_ends`, `_trigger`, `_lost`, `area_ledger_frac.extension`)
  in `present` and in every `series` entry, the explicit `params` the old run passed, and the `ext_*` / `extension`
  entries of `params_effective`. Every other `params_effective` value is equal. So the PARAMS defaults are the
  working baseline exactly.
- Pushed.

**Phase C** in `Infinite_Diffusion`: the uncommitted "Moved" banners were dropped, and the 299 files in
`migration/phase-c-remove.txt` were removed (`1e49bf2`). The README now points to this repository and the tag, and
`.gitignore` ignores personal `.cmd` files. A small regional `generate` wrote its GeoTIFF, and the GUI opened (user).
Pushed.

**The generator's tests had not run.** The before/after comparison saw "4 errors" both times. Four test files import
code that is not on `main`: `planet_diffusion.continents`, `.overview`, `.terrain_metrics` and `coarse.pool_coarse`,
from the `Solve_ScaleIssue_Attempt` branch work and the learned continent model of early October. Collection errors
stop pytest, so no test ran. Phase 0 had committed those test files to `main` because `tests/` held untracked files
from both branches. `Old Repo Tests.cmd` removed them from `main`, along with the scripts needing the same modules
(`measure_terrain_scale`, `verify_overview`, `verify_overview_region`, `compare_continent_controls`) and
`docs/continent-improvements.md`. It also dropped `--keep-checkpoints`, a branch-only flag (`main` always keeps
checkpoints), from `test_cube_workflow` and `test_draft`. Result: 152 passed (`2889f3f`, pushed). All removed files
stay in the tag. `Infinite_Diffusion` now tracks 63 files: `planet_diffusion/`, 17 tests, 4 scripts, 4 tools, 7 docs.

**Found on the way, not caused by the migration:** the source of `planet_diffusion/overview.py` and
`terrain_metrics.py` is in no branch, tag or folder. Only `__pycache__/*.cpython-312.pyc` from 2–4 Oct remain;
`continents` is on the branch. Copies of the three `.pyc` files are in `Infinite_Diffusion\outputs\recovered-pyc\`.

**State now:** both repositories are clean and pushed. Planet work happens only in `Planetary_Diffusion`. Model work
resumes with the rest of B (STATE section 5).
