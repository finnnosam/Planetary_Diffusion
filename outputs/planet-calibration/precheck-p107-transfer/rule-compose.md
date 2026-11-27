#### Question and reading rule for the composition-map probe (written before its runs, pass 107)

Finding that prompts it (paired traces from the same warm state, cloud seed 1028, 20 steps): the three return
transfers that treat the unread plate nodes ('conservative', 'matched', 'state') leave the end surface no better than
'read'. So the unread nodes are not the source. `returned` (the plate rasters read back at the same rotation right after
the transfer) already gains +87 holes on 'read': two world nodes read one plate node, the plate node can store one
value, and the nearest-node composition reads that node twice again. The map between two rotated copies of the same
grid is not one-to-one in either direction.

Question: does the non-bijective nearest-node read of the composition (with the return transfer that follows it) make
the per-step excess one-cell holes and crossings?

Arm traced: `compose_map` 'matched' (built for this probe: after ownership is decided by the nearest-node claims as
now, a world node that shares its plate node with a nearer world node of the same plate reads instead the nearest
unread live node of that plate among the shared node's eight neighbours; the transfer 'read' then returns each world
node's increment to the node it read). Against `default`, same warm state, 20 steps.

Rule:
- M1 (the build): matched's `returned` holes exceed its `transferred` by less than half of default's excess (+87).
  If not, the map is not one-to-one as designed: fix the build before reading M2.
- M2 (the premise): matched's mean `end` holes fall by >= 26 (25% of default's end 416 minus the smoothed start 313)
  and crossings by >= 2.5 per Mkm^2 (25% of the composition's per-step jump, ~10): the map is a source; it gets a
  full-length cloud pre-check and joins the sweep. Otherwise the source is the nearest-node resampling of node-scale
  structure at a new rotation whatever the map, or erosion's own one-cell features (eroded makes ~415 holes per step):
  the next question, not built this pass.
- The low flat interior and far low land are reported beside M2, not ruled on.
- One seed, a mechanism read; nothing adopted from it. The conservative return transfer (journal 143's arm) is swept
  regardless (built; the sweep decides).
