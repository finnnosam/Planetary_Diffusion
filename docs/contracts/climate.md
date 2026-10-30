# Contract (draft): climate (S3, S5)

**Nothing here is built.** This draft fixes the interfaces early so the plate and landscape stages can be made ready
for climate without waiting for it. It follows S3 and S5 of [`planet-ideal-flow.md`](../planet-ideal-flow.md); the
package is `planet_climate` (empty for now). Each item will be confirmed or changed when the stage is designed, with
its test set.

## S3 paleoclimate and sea level (along the plate history)

**Reads**, per epoch (every 1–5 Myr), from the [plate snapshots](plate-snapshots.md): `elevation_m`, `continental`,
`sea_level_m`, and the land / sea mask they imply; the plate run's `time_myr`.

**Writes** `outputs/planet-climate/<run>/paleo/t<time>.npz` on a coarse cube (n ≈ 54, ~185 km; final value set when
built), with the same time naming as the snapshots:

| Field | Unit | Consumer |
|---|---|---|
| `temperature_mean_c`, `temperature_range_c` (seasonal) | °C | S4 (frost, glaciers) |
| `precipitation_m_yr`, `runoff_m_yr`, `aridity` (P / PET) | m/yr, – | S4 (discharge, water balance) |
| `ice_thickness_m` | m | S4 (glacial erosion, ice load) |
| `sea_level_offset_m` (scalar: basin volume + ice) | m | S4 (shorelines, shelves) |

The landscape reads these by interpolating in time between epochs and bilinearly in space onto its own grid.

## S5 present climate and hydrology (after the landscape)

**Reads** the landscape's [`present.npz`](handoff-23km.md) (elevation, ocean, drainage) and latitude.

**Writes** `outputs/planet-climate/<run>/present.npz` on n = 434 (23 km; 11.53 km for the hydrology fields):

| Field | Unit | Consumer |
|---|---|---|
| `bio1` annual mean temperature, `bio4` temperature seasonality, `bio12` annual precipitation, `bio15` precipitation seasonality | WorldClim units | S6 bundle, S8 |
| monthly temperature and precipitation; PET; aridity index | °C, mm, – | S6, S9 |
| snowline / glacier equilibrium-line altitude; permafrost mask | m, – | S6, S8 |
| prevailing wind (u, v) | m/s | S6, S8 |
| discharge and lake levels from the water balance | m³/s, m | S6, S7, S9 |
| Köppen class, biome | class ids | S8, S9 |

**Earth side** (same operator rule): the S5 downscaler is trained on WorldClim 2.1, so every input it reads must be
computable on Earth from ETOPO by the same code.

**Checks** (from the ideal flow, bands to be fixed before the first run): latitudinal temperature gradient; desert-belt
latitudes (Earth ~20–30°); glaciated area through time; Köppen class shares (Beck et al. 2018); dryland share;
rain-shadow contrast across ranges.
