# Planetary Diffusion

A whole-planet terrain generator built from mechanism at large scales and learning at small ones:

```
plates and crust history (46 km, 250 Myr)  ->  landscape evolution (11.53 km, 30 Myr)
    ->  23 km hand-off  ->  learned meso model (2.88 km)  ->  Terrain Diffusion decoder (90 m, later)
```

Every stage is compared against Earth on the same grid (ETOPO, CRUST1.0, Bird 2003, GlobSed, GSRM, OCTOPUS,
Seton et al. 2020), with pass bands fixed before each run.

This repository was split on 29 Oct 2026 from
[Planetary-Infinite-Diffusion](https://github.com/finnnosam/Planetary-Infinite-Diffusion), which keeps the original
generator (Terrain Diffusion sampled on a sphere, GeoTIFF output). The history of the moved files came along; the
state before the split is that repository's tag `pre-split-2026-10-29`.

## Status

Research code, under active development; this README will become an install-and-use guide once the chain is
stable.

## Layout

One installable package per stage, under `packages/` (import names in brackets):

| Package | Stage |
|---|---|
| `planet-grid` (`planet_grid`) | cube-sphere grid, sphere graph and flexure, field helpers, constants, data/output paths |
| `planet-tectonics` (`planet_tectonics`) | S2 plates and crust history (46 km, 250 Myr) |
| `planet-landscape` (`planet_landscape`) | S4 landscape evolution (11.53 km) and the plate-to-landscape coupling |
| `planet-meso` (`planet_meso`) | S7 learned meso refinement (23 km to 2.88 km; needs PyTorch) |
| `planet-earth` (`planet_earth`) | Earth reference datasets on the same grids, and the data manifest |
| `planet-reviews` (`planet_reviews`) | reviews of the planet against Earth |
| `planet-climate` (`planet_climate`) | S3 / S5 climate (not built yet) |

`scripts/` holds the command-line runs and sweeps (`python -m scripts.simulate_plates`, `scripts.sweep_plates`,
`scripts.run_planet_landscape`, `scripts.planet_meso_check`, ...); `tests/` the tests.

- Design of every stage, its outputs and the statistics it is checked against: [`docs/planet-ideal-flow.md`](docs/planet-ideal-flow.md)
- What each stage writes and the next one reads: [`docs/contracts/`](docs/contracts/README.md)
- Where each parameter value comes from (measured, physical, assumed or fitted): [`docs/planet-parameter-audit.md`](docs/planet-parameter-audit.md)
- Background survey of planet-scale terrain methods: [`docs/research/planet-scale/README.md`](docs/research/planet-scale/README.md)

## Setup (Windows, Python 3.12)

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-dev.txt   # every package in editable mode, plus pytest
.venv\Scripts\python.exe -m pytest -q tests
```

The meso stage needs PyTorch (install the build for your GPU from pytorch.org first); the other stages need only
numpy, numba, ducc0 and matplotlib (rasterio for the Earth readers).

**Data.** The Earth datasets are not included. `python -m planet_earth.manifest` lists each one with its source and
the files the code reads, and checks the files you have against recorded checksums. They go in `data/`; results go
in `outputs/`. Neither folder is tracked. Set `PLANET_DATA` or `PLANET_OUTPUTS` to keep them elsewhere.
