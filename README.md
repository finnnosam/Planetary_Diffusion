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
stable. The layout is still the one inherited from the split (`planet_layer/`, `toy_meso/`, part of
`planet_training/`, `scripts/`); one package per stage is planned.

- Design of every stage, its outputs and the statistics it is checked against: [`docs/planet-ideal-flow.md`](docs/planet-ideal-flow.md)
- Where each parameter value comes from (measured, physical, assumed or fitted): [`docs/planet-parameter-audit.md`](docs/planet-parameter-audit.md)
- Background survey of planet-scale terrain methods: [`docs/research/planet-scale/README.md`](docs/research/planet-scale/README.md)

## Setup (Windows, Python 3.12)

For now the cube geometry comes from the original generator, so it is installed alongside:

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install numpy numba ducc0 matplotlib pytest
.venv\Scripts\python.exe -m pip install -e ..\Infinite_Diffusion --no-deps   # the original generator, checked out next to this one
.venv\Scripts\python.exe -m pytest -q tests
```

The meso stage also needs PyTorch. The Earth datasets are not included: they go in `data/` (ETOPO 2022, CRUST1.0,
Bird 2003 plate boundaries, GlobSed, GSRM v2.1, OCTOPUS, Seton et al. 2020; a download list with checksums is
planned). Results go in `outputs/`. Neither folder is tracked.
