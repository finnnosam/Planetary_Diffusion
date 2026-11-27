# Pass 108: remap_rule 'tent' — expectation on the trace (written before its trace runs)

Arm O2 of rule-floor.md (md5 842142d2…). The floor trace (trace-default-floor.jsonl, cloud seed 1028) read:
J_ret_next +172 holes / +9.5 crossings, J_world_next +13 / −1.4, J_ret_same +53 / +1.2; smoothing crust + sediment +
thermal age at `end` removes 57 / 46 / 61% of the excess (H / C / L). R4: the plate rasters' departure from the world
makes the jump.

Trace: the same warm state (−20 Myr), 20 traced steps with `--param remap_rule='"tent"' --floor`, against
trace-default-floor.jsonl. Expected (means over the 20 steps; default in brackets):

- composed_olddyn holes fall toward world_next's level: ~380–480 (573); crossings ~52–58 (63.5).
- rough_T's jump at the composition (+39 m, 601 → 640) falls below +10 m.
- end: holes ~300–380 (416), crossings ~46–51 (53.7), low flat interior ~5.5–7.5% of land (5.12), the gain growing
  over the 20 steps as the plate rasters lose the node-scale roughness the old transfer left.
- J_ret_next (the floor reading) under the arm: below half the default's (< +86 holes, < +4.8 crossings).
- If end's holes fall by less than 40 or crossings by less than 2 per Mkm² the arm does not do what R4 credited it with,
  and that is reported as such (the sweep still runs if the integrity pre-check passes: build and sweep, don't re-aim).
