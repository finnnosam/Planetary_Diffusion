# Pass 112: the frozen prediction for C's terrain claim (written 27 Nov before the Transfer Sweep started)

Change: C = `remap_rule` 'overlap' + `rift_area` 'plates' **with the material-aware return of pass 112** (journal 153)
against the default A = `remap_rule` 'nearest' + `rift_area` 'world'. Confirmation seeds 2001–2006, 250 Myr, saved
every 5 Myr, threads 2. Neither arm has been run on these seeds before. The development seeds 1001–1006 run in the
same sweep file (first call) and are read first; this prediction is not changed after they are read.

Read by `planet_reviews.geography` (the confirmation sweep's `reviews/geography.txt` / `.json`; n = 217, Earth ETOPO
without the ice regions): the six-seed mean of the paired differences C − A over the window −40 … 0 Myr (snapshots every
5 Myr); the endpoint 0 Myr reported beside it, not deciding.

Source of the ranges: the Rift Area Sweep's development seeds (sweep-20261126-161947, C before the repair; window paired
means, sd between seeds): holes −145 (37.7) of 358; low flat interior +3.02 (2.09) points; low land among far land
+5.47 (1.85); crossings −10.38 (6.06) of 53.35; land +0.24 (1.43); water 0–200 m deep −2.61 (0.59) points; shore step
+62 (56) m; low land among land far from the main water body +4.70 (4.09). A new six-seed mean differs from a
development one by about √2 × sd / √6 at one sd; the ranges are the development mean ± about 2.5 of that, rounded. The
repair changes where a return's loss is taken (arc material, not basement) and removes the uniform per-plate
correction; it is not expected to change these measures beyond that spread, which is why the old C's numbers are used.

| measure | predicted change, C − A, window, six-seed mean | deciding |
|---|---|---|
| one-cell water holes | −25 … −55% of A's | yes |
| low flat four-ring interior, % of land | +1.0 … +5.0 points | yes |
| low land among far land (≥ 300 km from any water), % | +3.0 … +8.0 points | yes |
| land–water crossings per Mkm² | −7 … −32% of A's | yes |
| land fraction, % of the globe | within −1.5 … +2.0 points (no claim on the sign) | yes |
| water 0–200 m deep, % of water (`shelf200`; C's cost) | no worse than −3.5 points | yes |
| shore step, median, m (`coast_step`; C's cost) | no more than +145 m | yes |
| low land among land far from the main water body (`low_far300_sea`) | positive | reported only |
| dry interior (four rings), far land (≥ 300 km), % of land | positive | reported only |
| land below 200 m, islands | no claim | reported only |

The prediction holds if all seven deciding lines fall inside their ranges. A failed line is reported, the terrain claim
is not made, and the failed line is the next thing to work on.

Beside it, scored, not deciding: the continents on the datum (`continent`, `level_split`, `three_factor`, both arms),
the hypsometry, and the vetoes V1–V5 on both seed sets as a check that C's integrity holds on new seeds; the oceanic
basement (`ocean_thin.py`: no oceanic node below 6,999 m in any snapshot of C, the invariant the repair restores).
