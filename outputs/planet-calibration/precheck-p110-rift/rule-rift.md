# Pass 110: the rift probe's question and reading rule (written before `rift_probe.py` runs)

The Remap Sweep (`sweep-20261126-141105`) failed V3 on C (`remap_rule` 'overlap'): the update part of the crust budget's
residual, whole run, is −0.36 km³/yr on C against −0.08 on the control (six-seed means; update / crystalline erosion
−0.123 … −0.177 against the control's −0.020 … −0.075, ±2 sd −0.083 … +0.003). Shift: −0.28 km³/yr.

Question: do the new plate nodes over continental rift gaps receive less crust than the world's pure shear takes from
the stretching zone (the zone strains by the gap's world cell area; the plates switch on the nodes whose forward nearest
node is the gap), by an amount of the order of that shift?

Probe: `rift_probe.py`, cloud seed 1028's warm state at −20 Myr (rebuilt this pass with the device's code), 10 steps
under 'overlap' and 10 under 'nearest' (same warm state; the hooks change nothing), and 10 under 'overlap' with
`rift_area` 'plates' (the arm built this pass).

Reading rule (per Myr, means over the 10 steps; 1 km³/yr = 10⁶ km³ per step):

- **R1 (the rift is the leak):** under 'overlap', shortfall (placed_rift − rift_world) ≤ −0.15 km³/yr and of the same
  size under 'nearest' (the rift's placement does not depend on the remap): the arm `rift_area` 'plates' is the
  next sweep's candidate, together with 'overlap'.
- **R2 (part of it):** −0.15 < shortfall ≤ −0.05: the arm is swept as built; the rest of the update part is read next
  (v3_stage's stages under C) and named in the journal.
- **R3 (not the rift):** shortfall > −0.05 km³/yr: the arm is still swept (built; the sweep decides), and the journal
  names the next term to read from the update part's gross terms (new oceanic-to-continental nodes, collision transfers,
  back-arc strips, arcs, holes).
- Under 'plates' the shortfall must be 0 to rounding (|·| < 10⁻⁶ of rift_world): a design check, not a result.
- other_gap (flow and erosion at gap nodes against the opened area) is reported beside it; it is not part of the rule.
