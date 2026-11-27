# Pass 111 probe: why the continents sit lower on the datum under the remap (rule written before the probe ran)

Seen in the Rift Area Sweep (sweep-20261126-161947, 1001–1006, 0 Myr; `level_split`): under C (overlap, plates) the
continents' lithosphere is thicker in every crust bin (mean L 30–35 km crust 126 → 151 km, Earth 107; 40–45 km
153 → 167, Earth 158), and the lithosphere term C carries the datum level's fall (Earth-weighted gap −29 → −230 m);
the Airy, density and sediment terms move little. L follows the thermal (lithosphere) age.

Question: is C's older-looking continental lithosphere

- H1, the remap's reading: the composition reads the thermal age by overlap weights (nine samples per world cell),
  so a narrow young zone (a reset by collision, arc or new crust) is shown mixed with its older neighbours, while the
  reset goes back to the plate nodes by 'read'; or
- H2, a different history: fewer or smaller resets over the run (collisions, arcs, rifts) under C.

Probe (`lith_age.py`): for both arms, every saved snapshot, six seeds: the area-weighted mean lithosphere age of
continental nodes, and the area share of continental nodes younger than 200 Myr; paired C − A by seed.

Reading rule:
- R1 (H1): the paired gap in mean age at the first saved snapshot (−245 Myr) is at least half its size at 0 Myr, with
  the same sign on at least five of six seeds.
- R2 (H2): the first-snapshot gap is less than a quarter of the 0 Myr gap.
- Between: both; reported as such. Either way the result names the next arm: R1 → the thermal age read by nearest
  plate node under 'overlap' (as the crust type and orogeny age already are); R2 → trace the resets' area per Myr by
  kind on both arms.
