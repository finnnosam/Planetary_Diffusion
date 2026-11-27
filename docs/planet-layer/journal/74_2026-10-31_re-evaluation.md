### Re-evaluation (31 Oct, after the fifty-second pass, at the user's request): the method is sound, but the work has drifted three steps from the terrain's largest failures into an open geodynamics question; proposed: park D's amplitude question, replace the fitted dynamic-topography field with Earth's measured spectrum, look at the whole chain's output on today's baseline, then return to the continental crust (too thick, on too little area); awaiting the user

**Why now.** The user asked for a whole-section double check ("in case we are getting lost in something we don't need
to"). The last direction change was pass 34 (journal 50, 26 Oct), which put B (continental budget and ocean basin) first.

**1. The goal, restated.** A whole-planet terrain generator: mechanism above ~100 km, learning below, judged against
Earth on the same grid. What the end user sees is elevation: hypsometry, mountain belts and plateaus, coasts and
shelves, the ocean floor, then the meso's 2.88 km texture. A pass is worth what it moves toward that.

**2. Where the effort went.**

| dates | passes | work | default path changed? |
|---|---|---|---|
| 17–18 Oct | steps 1–8 | the chain built and run end to end (plates, landscape, hand-off, meso) | built |
| 21–26 Oct | 1–34 | plate-layer calibration: sediment sink, shelves, breakup, Earth erosion law, margins, local isostasy, sediment layer, rift rule, Earth start crust, ring typing; extension E–E4 (5 passes) rejected | yes, many |
| 26–29 Oct | 35–41 | B, continental budget: B-a (compaction) and B-b (plate pieces) adopted; arcs, polarity, docking not adopted | yes (B-a, B-b) |
| 29–30 Oct | — | repository split and packages | no model change |
| 30–31 Oct | 42–44 | the intra-oceanic trench (T2): slab initiation and the torque law built, T2 "set aside as a direct lever" | no |
| 31 Oct | 45–52 | D: slab-driven mantle flow; then why it is 2.3–2.9× too strong on Earth (solver, geoid kernels, four viscosity profiles, SLB2011 slab thermodynamics, composition, the observation's sampling) | no (D off) |

**3. Signs of drift.**
- **Thirteen passes (40–52) without a change to the default path; the scorecard has not moved since 27 Oct.**
- **A chain of "the deeper cause is ...", each step reasonable alone:** too little continent and land too high (pass 34)
  → the continental budget (B) → the intra-oceanic trench excess (passes 36–44) → plates that move independently
  (pass 44) → a mantle flow to make them coherent (pass 45) → why that flow is too strong *on Earth* (passes 47–52).
  The work is now three steps from the failure it set out to fix.
- **The question passes 47–52 pursue is a known open problem.** Hoggard et al. (2016, Nature Geoscience) found
  observed long-wavelength dynamic topography about ±500 m against the ±2 km instantaneous flow models predict (~4×).
  Passes 47–52 reproduced it independently and carefully (E3 2.3–2.9×), and ruled out what the field has also ruled
  out (viscosity profile, slab buoyancy, sampling). The proposed next step (a 3-D lateral-viscosity flow) is
  research-frontier work, a large build, with no assurance of closing the gap.
- **The payoff for the terrain is small.** Earth's dynamic topography is ~0.22 km rms, 86% at l 1–3 (wavelengths over
  ~6,000 km). It tilts and floods broad lowlands; it is not what makes the planet look wrong. The visible defect
  journal 55 found (the 500 m, degree-30, plate-blind blobs: 20–30% of the ocean's high-pass variance) can be fixed
  by measurement directly: the audit's own step 5 ("amplitude and spectrum from one published global estimate, used
  as is") and journal 55's cheap control were never done.
- **Plate coherence (C(d), slow boundaries) is a real Earth difference, but a kinematic statistic.** It shows in the
  terrain only through what the plates build, which the crust measures already score.
- **The whole chain's output has not been looked at since 23 Oct.** The last landscape and meso runs
  (`outputs/planet-landscape/`, `outputs/planet-meso-check/`) are on the 22 Oct baseline, before the Earth erosion law,
  margins, isostasy, sediment layer, rift rule, start crust, typing, B-a and B-b. Since then the order has been set
  by scorecard lines alone, not by what the terrain gets wrong.
- **Decisions lean more on what was reported after the rule.** Pass 52's outcome (b) was resolved by a band split read
  after the score. That was honest and documented, but rules with three branches plus long "reported" lists give
  less protection than the simple rules of passes 20–34.

**4. What is right and stays.** The realism rule; values from Earth or mechanism; rules fixed before runs; six seeds;
Earth on the same grid; paired probes for small effects; step-by-step identity checks; the journal / decisions /
STATE discipline. D's code is verified (solver, geoid kernels, 'thermo') and stays built, off. Its Earth findings are
real and documented (journals 66–73).

**5. The largest failures, by what they do to the terrain** (scorecard, STATE section 4):
- land p25 / p50 / p90 858 / 1,444 / 3,317 m against Earth's 166 / 384 / 1,461: **land stands ~1 km too high**;
- land above 2 km 0.34–0.45 of land (Earth 0.05), and it is plateau (slope 2.0–2.5 m/km against 11.1);
- crust > 45 km 0.262 (Earth 0.086), > 55 km 0.075 (0.018); crust < 30 km 0.089 (0.190);
- continental area 0.321 (0.384), volume −9.9% over 250 Myr (≥ −3%), drowned continental crust 12.8% (~29%);
- sea level ~600 m below the datum.

These are mostly one problem: **continental crust that is too concentrated (too thick on too little area).** Pass 34
measured the split. Sea level explains about a quarter of the high land, and it is itself partly the lost area (more
deep ocean). The rest is the crust distribution. In STATE's order this is problems 3 and 4 (C, E), behind B and D.

**6. Proposed order (to agree with the user).**
1. **Park D's amplitude question.** No lateral-viscosity build. STATE records the finding as the literature's known
   discrepancy, reproduced; D stays built, off.
2. **Dynamic topography from Earth's measured spectrum** (audit step 5; journal 55 item 2): Earth's degree spectrum
   as published (Valentine & Davies 2020; Holdt et al. 2022 as a check), constant-variance blending between knots so
   it does not pulse. It replaces two fitted values (`dyn_std_m` 500, `dyn_slope` 1) with a measurement and removes
   the blobs. It stays plate-blind at l 1–3, a known, small, documented limitation. A small build, run as one arm of
   the next sweep. Later option: D's Earth-shaped pattern with a gain calibrated on Earth (E3 against Holdt's data,
   as the erosion law was calibrated on OCTOPUS). That changes the "no scale factor" rule, so it is the user's call.
   Not proposed now.
3. **Freeze and look** (STATE problem 7, moved up): one seed through landscape and meso on today's baseline, with
   hypsometry, maps at 65,536 colour levels, and the meso's flat-window check. This is to see what the terrain gets
   wrong now, and to set the order from it. Routine runs on the user's machine.
4. **The continental crust distribution (E with C).** Start with a measurement pass, no model change: where the
   > 45 km crust is and how it got there (collisions, the start, age of thickening), and where area is lost (ledgers
   by event). Then a mechanism not yet tried. E–E4 tried extension of thick crust and failed because thinning does not
   finish. Not tried: loss of the thickened lower crust (eclogitisation and foundering below ~50–60 km; there is
   no such term in the code); area lost at collisions (Scholl & von Huene's 0.7 km³/yr). Tests from pass 34's bands.
5. **B's remaining items** (T2, arc growth, sediment subduction) after that. The slab-flux excess mattered mostly for D.

**7. A working-rule suggestion** (the user's decision): a step-back check whenever five passes go by without a
default-path change, or every ten passes, against section 5's list. It asks whether the current thread still serves
the largest terrain failures.

Sources: Hoggard, M. J., White, N. & Al-Attar, D. (2016), "Global dynamic topography observations reveal limited
influence of large-scale mantle flow", Nature Geosci. 9, 456 (repository handle 1810/264180); Flament, N. et al.
(2013), "A review of observations and models of dynamic topography", Lithosphere 5, 189; Valentine, A. P. & Davies,
D. R. (2020), G3 21; Holdt et al. (2022), JGR 127; Scholl, D. W. & von Huene, R. (2009), Geol. Soc. Spec. Publ. 318.
