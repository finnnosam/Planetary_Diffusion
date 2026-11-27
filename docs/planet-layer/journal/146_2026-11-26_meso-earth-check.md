### Meso on Earth input (26 Nov, between passes 108 and 109; no runs): the current meso falls short of Earth's fine detail on hills and mountains even on Earth's own input; the meso is to be redone with more inputs (the user), so these become requirements for the new one

Asked whether an LLM could hand-edit the meso's weights to a good effect, and whether the meso had ever been tested on
Earth's own coarse input. Read from existing outputs only; nothing was run or changed in the model.

#### 1. What exists

- **October test** (`scripts/diagnose_meso_plains.py`, 15–16 Oct, before the repo split; the Project doc
  `claude/tectonic-prototype.md`): held-out Earth patches in three classes by 23 km relief (flat < 300 m, plain/upland
  300–800 m, hills 800–2,000 m), no mountain class. The adopted meso (`outputs/meso/detail-ft/step10000-passed.pt`) on
  real input, river incision / 5.8–11.5 km power, generated ÷ real: flat 1.12 / 1.15, plain/upland 1.10 / 1.69, hills
  0.96 / 0.94: passes the ≥ 0.7 rule. With blurred input and the true sub-cell detail map given: 1.16 / 1.58, 1.18 /
  2.03, 0.97 / 1.08; the planet supplies no detail map ("detail unknown").
- **Every planet meso check since 18 Oct** (`scripts/planet_meso_check.py`, column "meso on Earth input"): the meso
  run on the held-out split of `outputs/meso/dataset.npz` (boxes incl. Scandinavia, Tibet/Himalaya, Patagonia),
  windows classed by core p95–p5 relief, four classes. Identical in every run (same windows, same seeds;
  `results.json`, `verdict.<class>.earth_meso_*`):

| class | meso on Earth input: incision / 5.8–11.5 km power | meso on the planet, Chain run (age bug) / the fixed replay, band |
|---|---|---|
| flat | 1.06 / 1.28 | 1.72 / 2.79; 2.15 |
| plain/upland | 0.95 / 0.71 | 1.09 / 0.60 |
| hills | 0.83 / **0.46** | 0.72 / 0.35; 0.27 |
| mountains | 0.74 / **0.49** | 0.59 / 0.26; 0.34 |
| plains-lowlands | 0.85 / 0.74 | 1.46 / 1.02; 1.01 |
| coast | 0.96 / 0.80 | 1.78 / 0.87 |

#### 2. Reading

- On rough ground the meso falls short by itself: on Earth's own input it gives about half of Earth's 5.8–11.5 km power
  on hills and mountains (below its own 0.7 rule). About half of the planet's mountain-detail gap is the meso's; the
  rest is upstream (the landscape's missing high relief, journal 144 §1).
- The two Earth tests disagree on hills (fine-band power 0.94 in October, 0.46 in the planet check). They use
  different windows and class by different relief measures; not reconciled.
- Hand-editing weights is not a useful tool here. Measured edits (checkpoint blending, activation steering, a targeted
  fine-tune) could be scored honestly against Earth, but the meso is to be redone with more inputs (the user, 26 Nov),
  so no work goes into the current model.

#### 3. Requirements for the new meso (added to `docs/planet-ideal-flow.md`, S7)

- One held-out Earth test set covering all four relief classes, mountains included (the planet check's high-relief
  boxes), fixed before training; the current meso's numbers above are the baseline it must beat.
- Each candidate input earns its place the way the detail map did in October: on held-out Earth, the gain with the true
  value given against the input hidden, plus S6's planet-vs-Earth domain-gap check. An input the planet can only supply
  by fitting adds nothing at generation time.
- Which planet-side inputs are trustworthy (orogeny age, crust thickness, sediment cover, lithospheric support,
  rivers, climate) depends on §0's continental geography; the inputs are decided after it settles. Meso detail stays
  item 8 of STATE §5.

No code changed. Docs: this entry, `journal/README.md`, `decisions.md` (one row), STATE (§4a meso line, §5 item 8),
`docs/planet-ideal-flow.md` (S7 stats; published), `next-chat-prompt.md` (journal counter).
