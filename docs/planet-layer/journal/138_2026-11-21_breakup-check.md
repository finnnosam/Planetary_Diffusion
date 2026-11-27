### Breakup check (21 Nov, after pass 102; the user saw seed 1001's supercontinent stay whole under 'earth_corr'): breakup is weaker under the adopted field on 5 of 6 seeds; seed 1001's supercontinent never splits; the field enters neither the forces nor the split rule, so the path is indirect and unknown; pass 103 reads it first

Probes: `outputs/planet-calibration/precheck-p102-share/plate_continent.py` (printout `plate_continent.txt`),
`breakup_pieces.py`. Sweep `sweep-20261121-185217` (journal 137).

#### 1. What the user saw, measured

`flow_placed-earth_corr` seed 1001: the start's supercontinent plate (plate 0) keeps 182 of the planet's 208 Mkm² of
continent to 0 Myr (label-free, crust + bulk sediment ≥ 20 km), and it is one connected piece. Its rifts shed only
19.2 (−161 Myr), 1.8, 1.3, 2.4, 1.4 and 3.3 Mkm²; in the control's seed 1001 plate 0 split off 88.1 Mkm² at −153 Myr
and kept splitting (24 rifts against 15). Its rifted production is 610 km/Myr (the arm's other seeds 1,773–2,303;
Earth 1,553).

The most continent held by one plate at 0 Myr, by seed:

| | 1001 | 1002 | 1003 | 1004 | 1005 | 1006 | mean |
|---|---|---|---|---|---|---|---|
| A `flow_placed-all` | 96 | 91 | 37 | 52 | 51 | 62 | 65 |
| B `flow_placed-earth_corr` (adopted) | 182 | 67 | 61 | 96 | 114 | 104 | 104 |
| pass 100's model (`sweep-20261121-133838`, the other trajectory) | 96 | 119 | 38 | 43 | 39 | 116 | 75 |
| pass 98's model (`dyn_rule` 'measured') | 74 | 46 | 55 | 92 | 61 | 77 | 67 |
| Earth (Bird 2003 plates on CRUST1.0) | | | | | | | 42 (then 35, 27) |

The largest plate at 0 Myr is larger under B on 6 of 6 seeds (126–195 against 72–131 Mkm²); B − A per seed in
continent on one plate +86 / −24 / +24 / +44 / +63 / +42 (mean +39, se ~16). The two trajectories of the same model
differ by 10 in the mean: B's +39 is beyond that range, though six seeds with an sd near 40 do not settle it.

Two readings that do not separate the arms: the continent's connected pieces on the graph (largest 150–190 of ~205
Mkm² on most runs, both arms) are Earth's kind too (Earth's continent on the same grid and rule: 169 of 188 Mkm² in
one piece, joined by shelves and arcs), so connectivity is not the measure; the plate is.

Already before B the planet held 1.5–1.8× Earth's continent on one plate (65–75 against 42): weak breakup is an old
failure (STATE §5 item 3, "the breakup's driver"; the passive edge 0.55×), and B made it worse.

#### 2. Where the field could enter

Read in the code: the dynamic field enters neither the boundary forces that rank the rift lines (`_cut_forces`: slab
pull, slab suction, ridge push, drags; no elevation term) nor the split rule (`_split_events_earth`: Earth's rifted
production × continental area per km of the candidate's margin) nor the plate motions (the flow's velocity only with
`basal_traction`, off). So the field reaches breakup only through the state: elevation → erosion, sediment, crust
thickness and the continental flag; which crust drowns; where subduction starts (B logs 'initiation at margin
skipped' 3× on seed 1001, A none); and plain divergence of trajectories. Two features of the split rule make a whole
supercontinent hard to split whatever the cause: the acceptance per step goes as 1 / L (a long, plate-wide cut is
accepted far less often than a sliver), and with `rift_line` 'force' the candidates are ranked by opening force, which
a plate ringed by trenches it overrides (slab suction fitted to 0) lacks except where an attached slab pulls.

#### 3. What follows

Adoption stands by the rule fixed before the runs (breakup is scored, not a veto; the setting is one default, cheap to
undo if a later reading shows the cost is the field's). But an unbroken supercontinent bears on the target itself (one
huge interior, its highs and its drowned rim), so **pass 103 first reads why B breaks up less**, before building the
field–continent arm: per seed, the time series of plate 0's area and of the candidate cuts' kind, opening force and
margin length (the chord / bisector / force counts in `stats.json` already differ: 171 / 262 / 890 against 264 / 158
/ 1,055 on seed 1001), and whether the state changes the field makes (continental flag, trench placement) change the
cuts. If the cause is the rule's 1 / L bias or the missing opening force, the breakup's driver moves up §0's order
(to the user, as a change of order).

Scorecard line added to STATE: continent on one plate at 0 Myr, mean of seeds.
