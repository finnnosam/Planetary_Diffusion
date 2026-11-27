# Pass 108: the floor's second trace — question and reading rule (written and committed before any trace runs)

Journal 144 §5. Cloud, seed 1028, the warm state rebuilt to −20 Myr by `trace_shore.py warm` (pass 107's script,
unchanged), then 20 traced steps of the default model with `trace_shore.py trace --mode default --floor` (this
folder's copy). The new readings are side readings on each step's `end` state (the trajectory is the untraced one:
`check-identity` with the readings on). One seed: a mechanism read, not a terrain claim.

## Question

Which material roughness, under which support, makes the end-of-step excess of one-cell holes, land–water crossings
and the missing low flat interior; and is the composition's jump made by the nearest-node resampling of the world's own
surface, or by the plate rasters' departure from it?

## Readings (each a full surface on the world grid, its own flexural or local-Airy solve, the sea re-solved and restored)

Q1, the factorial at `end`. Smoothing = conservative graph diffusion over the world's continental nodes (not this
step's new crust), 100 km as pass 106 (`smooth_material`'s length), on the world state at `end`:
- `F0` = `end`; `A0` = `end_local` (local Airy in place of flexure, as pass 105's side reading).
- `F_T`, `F_S`, `F_A`, `F_TS`, `F_TSA`: flexure, with the crystalline crust (T), the sediment (S), the thermal age
  `therm` (A: the C term's input), T and S together, and all three smoothed.
- `A_T`, `A_S`, `A_A`, `A_TS`, `A_TSA`: the same under local Airy.

Q2, the composition at the same and the next rotation (all with `end`'s dynamic field and sea):
- `ret_same`: the plate rasters read at this step's rotation with this step's ownership (after the events).
- `ret_next`: the plate rasters read at the next step's rotation (each plate's Q advanced by its current axis and speed,
  as `step` does), ownership by `_claims` + `_world_from` at that rotation (state saved and restored), t + dt.
- `world_next`: the world's own `end` state carried to the next rotation by the same nearest-node rule: a world node
  owned by plate k at the next rotation reads the `end` world node nearest Q_k Q'_kᵀ x (where its material is now) when
  that node is owned by k and not new crust, else the `ret_next` value (fallback share recorded).
- `world_next_local`: `world_next` under local Airy.
- Validity: `ret_next` must reproduce the next step's `composed_olddyn` (mean holes within ±10%, crossings within ±5%).

Measures (means over the 20 steps): holes H, crossings C (per Mkm²), low flat four-ring interior L (% of land); Earth on
the same grid (`geography.py`): H 13, C 19.14, L 14.24%. Excess e(M) = M(F0) − M(Earth) for H and C, and
M(Earth) − M(F0) for L. The share a reading D removes: s_D(M) = (M(F0) − M(D)) / e(M) (sign reversed for L).

## Reading rule

- R1 (support). Report s_A0 on H, C, L.
- R2 (material under flexure). A field X ∈ {T, S, A} *makes* part of the floor if s_FX ≥ 0.25 on at least two of H, C,
  L. If s_F_TSA ≥ 0.50 on two of three, the materials' node-scale roughness is the floor's main source. A smoothing that
  makes the surface rougher (s < −0.10 on two of three) marks a field that compensates another (reported).
- R3 (channel). For a named X, its share under local Airy measured from A0 on the same excess, s'_AX = (M(A0) − M(AX)) /
  e(M). If s_FX ≥ 2 s'_AX on the measures that named it, the channel is flexure's support of X's node-scale roughness.
- R4 (resampling). Jumps J_D(M) = M(D) − M(F0) for H and C. If J_world_next ≥ 0.5 J_ret_next on both, the resampling of
  the world's own surface makes most of the jump (journal 144's suspect a). If J_world_next < 0.25 J_ret_next on both,
  the plate rasters' departure from the world makes it. Otherwise both (reported with the numbers). J_ret_same is the
  departure at the same rotation (shared reads).

## Which arm follows (fixed now)

- O1: A named by R2 (the thermal age): the C term's buoyancy as a load on the elastic plate (today it stands locally
  compensated at every node, the Te → 0 limit; the thermal anomaly lies in or under the plate). Built off, swept.
- O2: T or S named by R2 (alone or in TS), or only TSA jointly (≥ 0.50): the materials' node-scale roughness is
  numerical (the transfer and the resampling between rotated grids of one resolution). Arm: a partition-of-unity kernel
  remap between the world and the plate rasters (the composition gathers each world node from the plate nodes around
  Qᵀx with normalised tent weights; the increments go back by the same weights, conserved per plate), so no plate node
  is unread and no node-scale structure is copied a node off. R4 states what it should remove (both parts if the jump
  is the world's resampling; the departure if it is the rasters'). It accumulates state: full-length cloud pre-check
  (seed 1030) before its sweep.
- O1 and O2 together: both are built; the larger s_F (summed over H, C, L) is swept first.
- O3: none named and s_F_TSA < 0.50, but s_A0 ≥ 0.40 on two of three: the support matters through something not
  smoothed here (the ocean floor at the coast, the column water-load decision, sharp continental edges). The pass is an
  evidence checkpoint and names that split as the next trace.
- O4: none of the above: the floor is erosion's own (the landscape law at 46 km); evidence checkpoint, reported to the
  user with STATE §0 item 1's order in question.
