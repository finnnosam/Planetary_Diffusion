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

Research code, under active development. The layout is still the one inherited from the split
(`planet_layer/`, `toy_meso/`, part of `planet_training/`, `scripts/`); restructuring into one package per stage
is planned (`docs/planet-layer/STATE.md`).

- Current state, scorecard against Earth and next steps: [`docs/planet-layer/STATE.md`](docs/planet-layer/STATE.md)
- What was adopted or rejected, and why: [`docs/planet-layer/decisions.md`](docs/planet-layer/decisions.md)
- Code, scripts, data and commands: [`docs/project-map.md`](docs/project-map.md)
- Long-term design of all stages: [`docs/planet-ideal-flow.md`](docs/planet-ideal-flow.md)

## Setup (Windows, Python 3.12)

For now the cube geometry comes from the original generator, so it is installed alongside:

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install numpy numba ducc0 matplotlib pytest
.venv\Scripts\python.exe -m pip install -e ..\Infinite_Diffusion --no-deps   # the original generator, checked out next to this one
.venv\Scripts\python.exe -m pytest -q tests
```

The meso stage also needs PyTorch. The Earth datasets are not included: they go in `data/` (see
`docs/project-map.md` for each dataset and where it comes from). Results go in `outputs/`. Neither folder is
tracked.
