# Pass 111: the frozen prediction for the remap's terrain claim (written 26 Nov before the confirmation runs)

Change: `remap_rule` 'overlap' with `rift_area` 'plates' (arm C, PARAMS since pass 111 as a numerics fix) against the
old default `remap_rule` 'nearest' with `rift_area` 'world' (arm A). Confirmation seeds 2001–2006, 250 Myr, saved every
5 Myr, threads 2. Neither arm has been run on these seeds before.

Read by `planet_reviews.geography` (the sweep's `reviews/geography.txt`; n = 217, eight-neighbour graph, far land by
graph distance, Earth ETOPO without ice): the six-seed mean of the paired differences C − A over the window −40 … 0 Myr
(snapshots every 5 Myr); the endpoint 0 Myr reported beside it, not deciding.

Source of the ranges: the Rift Area Sweep's development seeds 1001–1006 (sweep-20261126-161947, the same pair of arms):
window paired means (sd between seeds) holes −145 (37.7) of a control 358 (−40%); low flat interior +3.02 (2.09)
points; low land among far land +5.47 (1.85); crossings −10.38 (6.06) of 53.35 (−19.5%); land +0.24 (1.43) points;
land below 200 m +1.33 (3.47); dry interior +7.06 (6.28); far land +7.20 (6.97). A new six-seed mean differs from the
development one by about √2 × sd / √6 at one sd; the ranges are the development mean ± about 2.5 of that, rounded.

| measure | predicted change, C − A, window, six-seed mean | deciding |
|---|---|---|
| one-cell water holes | −25 … −55% of A's | yes |
| low flat four-ring interior, % of land | +1.0 … +5.0 points | yes |
| low land among far land (≥ 300 km), % | +3.0 … +8.0 points | yes |
| land–water crossings per Mkm² | −7 … −32% of A's | yes |
| land fraction, % of the globe | within −1.5 … +2.0 points (no claim on the sign) | yes |
| dry interior (four rings), far land (≥ 300 km), % of land | positive | reported only |
| land below 200 m, islands | no claim | reported only |

The prediction holds if all five deciding lines fall inside their ranges. A failed line is reported, and the terrain
claim is not made; the change stands as the numerics fix it was adopted as (pass 111), and the failed line is the next
thing to work on.

Beside it, scored, not deciding: the continents on the datum (`continent`, `level_split`, `three_factor`, both arms),
the hypsometry (`hypsometry` reviews), the vetoes V1–V5 as a check that C's integrity holds on new seeds.
