# Pass 112: the rule, written before the repair is run (27 Nov 2026)

## Question

Does a material-aware return remove C's unbooked thinning of oceanic basement, without a clamp doing the work and
without touching the default path?

C = `remap_rule` 'overlap' + `rift_area` 'plates'. The repair changes only the 'tent' / 'overlap' return
(`composition._tent_return`):

1. Each world increment is split into its gain and its loss before it is mixed onto plate nodes, so a cut next to a fill
   cannot cancel and then be amplified (the pass 109 blow-up came from splitting *after* mixing).
2. Gains go back by the overlap weights, as before (each plate node the area mean of the gains over it).
3. A loss on a world cell is debited from its contributors in proportion to weight × available material, where
   available = the node's value above its floor (oceanic thickness: `ocean_thick_m`, the floor the forearc cut itself
   respects; continental thickness, sediment, rock: 0). A contributor with nothing above the floor takes no debit.
4. Each plate's gains and losses, per crust type, are closed to the world's by rescaling **only the nodes that took that
   sign** (no uniform correction over the plate, which moved untouched basement: journal 151 §3).
5. After rescaling, a loss larger than a node's available material is cut back to it, and the volume not removed is
   booked (`erosion_totals['return_shortfall_<field>_m3']`). Loss and gain that no live node can take are booked too.

## Pass bands and the rule that decides

- **R1 (unit test on a physical example):** a basement patch at 7 km beside an accreted arc; a known forearc cut on the
  world: every basement node unchanged (exactly), removed = requested on the plate within float32 rounding, shortfall 0;
  a continental cut in the same call is conserved on its own; repeated with a node no overlap reaches (takes the world
  node over it, debited only against its own material) and with an ownership change (an arc contributor no longer owned:
  its share goes to the other contributors with material, or is booked, never onto basement).
- **R2 (the trace, `trace.py`, the review's `remap_trace.py` adapted to the new return):** 10 steps on seeds 1001 and
  1002 and 250 steps on 1001: **zero** oceanic nodes that cross below 6,999 m in the return without the volume being
  booked; the run's oceanic deficit below 7 km ≤ booked oceanic shortfall-free losses (i.e. the deficit at 0 Myr is 0 or
  is fully booked).
- **R3 (the cap is not the mechanism):** the booked shortfall over the 250 Myr trace ≤ 0.1% of the loss volume returned.
  Above that, the cap is acting as a clamp: not accepted, traced further.
- **R4 (the suite with C as the default):** the two invariants that failed in pass 111 pass (oceanic minimum ≥ 6,999 m;
  oceanic deficit ≤ booked oceanic erosion); the three trajectory checks named in journal 150 §3 are pinned with notes
  if they still fail; the six option-setup tests are given 'nearest' explicitly.
- **R5 (default path):** the change is inside the 'tent' / 'overlap' branch only; the planet suite passes with the
  defaults and `precheck-p104-repeat/rep.py` (n = 48, 70 steps, threads 1 and 2) matches the committed default.

If R1–R5 pass: C (repaired) is swept against the control on 1001–1006 and 2001–2006 in one sweep, with a terrain
prediction frozen first (from this repaired C's development-seed trace window, template
`precheck-p111-confirm/prediction.md`, adding `shelf200` and `coast_step`). If R2 or R3 fails, the pass ends as an
evidence checkpoint naming what is left.
