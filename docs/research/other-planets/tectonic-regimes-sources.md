# Sources of "An Apple Pie From Scratch, Part Vb" (tectonic regimes): triage for other planet types

> **Priority: low (user, 7 Nov 2026).** Other planet types are a later goal, not the current work. This note
> supplements the Earth model and does not change its order. Use it for background, and for the Earth-relevant
> items (Mallard 2016 plate statistics, subduction-initiation papers) where they fit a pass already planned.
> Don't start other-planet work from it.

Written 7 Nov 2026 (chat, not a pass; nothing in STATE or decisions changes). Article:
https://worldbuildingpasta.blogspot.com/2020/01/an-apple-pie-from-scratch-part-vb.html (worldbuilding blog, updated
through 2025; 134 numbered references). Context: the user wants other surface types eventually (Martian, Venusian, ...).
This note says which of the article's sources help with that, and which help the Earth model now.

## What the sources are, in one paragraph

Part Vb is about the *interior regime* (stagnant, squishy, sluggish, episodic, mobile lid) and its link to climate
and life. Its sources are mostly geodynamics models and early-Earth geology, plus a block of habitability papers. It
cites **no topography datasets, no crater production or degradation laws, no hypsometry studies and no surface-process
(erosion, aeolian) papers**. So it answers "which regime is this planet in and why", not "what does its surface look
like". Of the model papers checked, **none report topography or elevation statistics**: they report plateness,
mobility, crust and lithosphere thickness, heat flow and surface age.

## Spot-checked (read in full or in part, 7 Nov)

| Ref | Paper | What it gives us | Kind |
|---|---|---|---|
| 86 | Mallard et al. 2016, *Nature* 535 | 3D StagYY, yield stress 100–250 MPa. Plate areas fall in two populations (≈7 large plates + a fractal small-plate tail); crossover ≈ 10^7.6 km² (≈ 4×10^7 km²) for Earth (Bird 2003) and the best model. Largest plate size follows slab spacing (convection degree 6 → 4 → 2 as yield stress rises). Triple junctions sit more on trenches than ridges: Earth 131 vs 71. Plate-like behaviour 100–350 MPa. | Model + Earth measure |
| 45 | Byrne, Foley et al. 2021, *JGR Planets* 126 | Brittle–ductile transition depth from Byerlee friction + olivine dislocation creep and Peierls plasticity, ~10^4 models over g 5–40 m/s², T_s 200–1000 K, T_m 1500–2000 K, plate age 5–300 Myr, strain rate 10^-13–10^-16 /s. Surface temperature is the strongest control (T_s 1000 K → ≤ 5 km); higher g → shallower; mantle temperature negligible; hydrostatic pore pressure deepens it 5–10 km at low T_s and g. Earth ≈ 25 km (100 km lithosphere), Mars ≈ 40–60 km. Thin brittle layers predict little high topography ("eggshell"). **MATLAB code and outputs on Zenodo** (10.5281/zenodo.5560138, 10.5281/zenodo.5560142). | Mechanism |
| 71 | Lyu, Ballmer et al. 2025, *Nat. Commun.* 16 | 2D half-annulus StagYY with melting, 10 Gyr runs. Six regimes, including a new **episodic-squishy lid** (proposed for Archean Earth and Venus). Mobile lid for yield stress ≲ 90 MPa at any thermal state; hot: episodic 130–180 MPa, stagnant ≳ 180; moderate: squishy ≳ 120; cool: sluggish ≳ 120. Authors: thresholds are not universal (depend on the viscosity profile). Intrusion:extrusion fixed at 9:1. Model data on Zenodo (10.5281/zenodo.14601956); StagYY not public. No topography. | Model |
| 34 | Lourenço et al. 2020, *G-cubed* 21 | 2D annulus StagYY. **Plutonic-squishy lid** when intrusion efficiency ≳ 70% and surface yield stress ≳ 50 MPa: small, strong, short-lived plates (boundaries live a few Myr), moving up to several cm/yr without subduction, driven by delamination and eclogite drips; thin lithosphere; crust 50–100 km; high heat flow. No topography. | Model |
| 43 | Tian, Tackley & Lourenço 2023, *Icarus* 399 | 2D, experimental plagioclase (An75) crust rheology + intrusive magmatism. Plagioclase crust → "deformable episodic lid" (distributed strain between overturns); olivine crust → episodic or stagnant. Best Venus match: plagioclase crust, low eruption efficiency (high eruption gives surface ages < 100 Myr). Abstract only. | Model |
| — | Rolf et al. 2018, *Icarus* 313 (fetched in place of ref 70, Uppalapati 2020, by the same group) | 3D StagYY Venus; geoid and dynamic-topography spectra to degree 16 compared with Venus gravity (SHGJ180U) and topography; episodic (yield 50 MPa) fits better than stagnant (misfit 0.26 vs 0.52); mean surface age 0.25–0.6 Gyr; overturns mobilise the surface for 150–200 Myr. The one checked model that tests long-wavelength topography. | Model vs Venus data |
| 37 | Ghail et al. 2024, *Space Sci. Rev.* 220 | Venus observational numbers (first ~15% of the review read): ~1,000 craters, ~10% modified; mean surface age 150–1,000 Ma; volcanic plains ~92%, tessera ~8%; ~500 coronae, 60–2,600 km, median ~200 km, ~62% along rifts or fracture belts; > 85,000 volcanic edifices (Hahn & Byrne 2022: ~32,500 < 5 km, 729 at 5–100 km, 118 > 100 km); hotspot rises 1,000–2,000 km wide, 1–2 km high; largest volcanoes ~9 km; Magellan altimetry and stereo too coarse for volcano slopes. | Observation (a Venus scorecard) |
| 35 | Byrne, Ghail et al. 2021, *PNAS* 118 (read in full from the user's PDF, 7 Nov) | **58 campi**: radar-smooth lowland plains ≥ ~50% enclosed by ridge belts (shortening, folds over thrusts) and/or groove belts (graben). Sizes ~100 to > 1,000 km; the largest is 1,700 × 1,100 km (≈ Alaska), probably near the largest block Venus's thin strong upper crust can hold together. 79% lie < 1 km above the 6,051 km reference radius; 57% sit on geoid lows. Interiors barely deformed (wrinkle ridges); margins carry a few % strain (example: 2.6–3%, 2.76–2.8 km extension over 106 km), with transpressive/transtensive shear, so blocks translated and rotated, some after the youngest plains. A weak lower crust at ~5–20 km (Earth's equivalent ~50 km) decouples them. Driver: stresses from gravity-inferred mantle flow (propagator-matrix flow kernels, radial viscosity from dry diabase crust and dry dunite mantle; nominal lithosphere 15 km, 18 mW/m², 10^-15 /s) exceed Byerlee failure in the upper few km at every campus. The authors link it to Lourenço 2020's plutonic-squishy lid. Mapping: Magellan SAR mosaics (~75 m/px) at 1:15,000,000, altimetry 4.6 km/px. Per-block areas, centroids, elevations and geoid values are in SI Table S1 (not in the main PDF). | Observation + mechanism |

## The rest, by use (not opened)

**For the Earth model now** (mechanism or scored measures):
- 97 Cawood & Hawkesworth 2019: continental crustal volume, thickness and area (compilation of measurements).
- 57 Korenaga 2013; 64 Gerya 2014: reviews of plate-tectonics initiation and Precambrian geodynamics.
- 58 Nikolaeva et al. 2010 (subduction initiation at passive margins); 59 Ueda et al. 2008 (plume-induced
  initiation); 53 Yuan et al. 2024; 65 O'Neill et al. 2017; 60 Ruiz 2011 (impact-triggered): subduction-initiation
  mechanisms.
- 63 Rolf & Tackley 2011 (continents focus stress in 3D convection); 84 Cathles et al. 2023 (asthenosphere's role);
  103 Lenardic et al. 2005 (continental insulation and ocean/continent area); 98 Höning et al. 2019 (bifurcation in
  continental growth); 99, 100 Huang & Li 2022, 2023 (supercontinent assembly and lifetimes); 85 Yoshida 2023.

**Regime mechanics** (what would move a generated planet into another regime):
- Heat pipe: 15 Moore & Webb 2013; 17 Moore et al. 2017; 61 Kankanamge & Moore 2016; 16 Reinhold & Schaefer 2025.
- Regime transitions and diagrams: 74 Weller & Kiefer 2020; 76 Foley 2018; 79 Lenardic 2018; 80 O'Neill et al. 2016;
  87 Xiang et al. 2025 (3D modes); 72 Lourenço et al. 2016; 36 Rozel et al. 2017; 38 Jain et al. 2019; 81 Rozel et al.
  2015 (ridges in a stable lithosphere); 82 Capitanio et al. 2019 (lid-and-plate); 78 Al Asad & Lau 2024 (sluggish lid).
- Planet size and composition: 47 Korenaga 2010 (super-Earths); 51 Skemer et al. 2025 (composition and shear
  localization); 50 Unterborn et al. 2017; 52 Noack et al. 2014; 46 Putirka 2024.
- Climate → tectonics: 48 Lenardic et al. 2008; 49 Landuyt & Bercovici 2009.
- Hemispheric: 89, 90 Meier et al. 2024, 2021; 88 Van Summeren et al. 2011; 91 Citron et al. 2018 (Mars dichotomy:
  degree-1 convection antipodal to an impact); 92 Zhong et al. 2000 (lunar mare asymmetry).
- Venus: 19 Gerya 2014 (plume-induced crustal convection → coronae/novae, 3D); 66 Davaille et al. 2017
  (plume-induced subduction at coronae); 41 Capitanio et al. 2024 (Ishtar Terra); 42 Montiel et al. 2025 (chasmata);
  40 Collinet et al. 2025; 70 Uppalapati et al. 2020 (lithospheric overturns); 106 Yin 2012 (Tharsis by slab rollback,
  a minority model).

**Observations** (possible test data for other bodies):
- Venus: 18 Smrekar et al. 2023 (lithosphere thickness and heat flow at coronae); 22 Strom et al. 1994 (crater
  distribution and global resurfacing); 25 Filiberto et al. 2020.
- Io: 12 McEwen et al. 2004; 13 Mura et al. 2025 (Juno lava lakes); 14 Schenk & Bulmer 1998 (mountains by thrusting).
- Moon: 33 Solomon & Head 1980 (mascon basins); 30, 31, 32 (young volcanism ages).
- Mars: 26 Horvath et al. 2021; 27 Kite & Conway 2024.
- Icy bodies: 127 Stern, Gerya & Tackley 2018 (stagnant-lid review across bodies); 128, 130, 131 Europa; 132 Barr &
  Hammond 2015; 133 Liu et al. 2016 (Titan structural map); 129 Collins et al. 2009.

**Not about terrain** (habitability, climate, oxygen, carbon cycle, Hadean geochemistry): 1–10, 23, 24, 28, 54, 56,
67, 68, 77, 83, 101, 102, 104, 105, 107–126.

## Quality flags

- **One code, one group.** Most regime models (Tackley, Lourenço, Rozel, Jain, Tian, Lyu–Ballmer, Rolf, Mallard)
  use StagYY, mostly in 2D annuli at Earth-scaled or lower Rayleigh numbers. Agreement between them is not
  independent confirmation; Lyu et al. say their thresholds depend on the setup.
- **Minority or contested views**, cited by the article as the edge of a range: 55 Hamilton 2011 (plate tectonics only
  since the Neoproterozoic, no deep plumes); 73 Piper 2018 (lid tectonics from palaeomagnetism); 69 Stern 2025
  (two plate and two single-lid episodes); 106 Yin 2012.
- **Not peer-reviewed papers:** 9 (EGU abstract), 134 (LPSC abstract); 29, 50 (arXiv preprints); 85, 93, 107, 129
  (book chapters); 20, 95 (no venue given). The Notes also link a creationist PDF (icr.org), cited as a joke, not
  as a source, and a YouTube talk on tectonic modes and climate (youtube.com/watch?v=ak-IKxgaZ5Q).

## What this means for the project

1. **Regime is a mechanism question; surface form is an observation question.** The model literature can give a
   generated planet its regime, crust and lithosphere thickness, and resurfacing age from physical inputs. What the
   surface looks like must be scored against each body's own data (as ETOPO is for Earth), and the article cites none.
   For each body, the next research step is the elevation data and the statistics the decisions would rest on:
   hypsometry, crater size–frequency, feature catalogues.
2. **Directly reusable mechanism:** Byrne et al. 2021 maps (g, T_s, plate age, strain rate) to brittle-layer thickness
   with published code; Lyu 2025 and Lourenço 2020 give regime boundaries in yield stress and intrusion efficiency.
   These are the inputs a "planet type" switch would turn.
3. **For the Earth model now:** Mallard 2016's plate-area crossover (≈ 4×10^7 km²) and triple junctions by boundary
   type (Earth 131 on trenches vs 71 on ridges) are Earth measurements not yet on the scorecard (it has the plate
   count, ≤ 52, and plates ≥ 1 Mkm²). Its mechanism (largest plate size set by slab spacing) is a cross-check on the
   plate code. Candidates for a scored line, not vetoes.
4. **Venus scorecard seed:** Ghail 2024's counts (coronae, volcanoes, tessera fraction, rise sizes, crater count),
   plus Byrne 2021 (PNAS)'s 58 belt-bounded lowland blocks (statistics below).

## Byrne 2021 (PNAS) SI Table S1: the 58 campi (extracted 7 Nov from the user's SI PDF)

Per block: centroid latitude and longitude, area (nearest 1,000 km²), centroid elevation (Magellan, relative to
6,051 km radius), geoid, and modelled von Mises stress at 5 km. Saved as
`docs/research/other-planets/venus_campi_tableS1.csv` (58 rows, ids 1–58 all present).

| Measure | Value |
|---|---|
| Total area | 16.3 Mkm², **3.5% of Venus** (they are a feature catalogue, not a tiling of the surface like Earth's plates) |
| Area | 14,000 to 1,716,000 km²; median 150,000; mean 280,000; largest 5 hold 42% of the total |
| Equivalent diameter (2√(A/π)) | 134 to 1,478 km; quartiles 291 / 436 / 610 km |
| Size distribution | Cumulative N(>A) ∝ A^−0.9 to A^−1.1 for A ≥ 10^5 km² (least squares, 36 blocks; maximum likelihood 0.92 ± 0.15). Flatter below ~10^5 km² (0.66 from 4×10^4), most likely because small blocks were missed at the 1:15,000,000 mapping scale |
| Centroid elevation | −556 to +2,137 m; median 263 m; 79% below 1 km (matches the paper); 24% below the reference radius |
| Geoid | 55% negative, 59% ≤ 0 (paper: 57%) |
| Stress at 5 km | 22–198 MPa, median 86 |
| Area vs elevation | Barely related (r = 0.19 between log area and elevation) |
| Hemispheres | 30 north, 28 south |

Use as a test: a Venus-regime run would be scored on the count of rigid low-lying blocks bounded by deformation
belts, their size distribution above ~10^5 km² (slope about −1), and their elevation (mostly within 1 km of the
mean radius). Caveats: one mapping team's definition (≥ 50% of the margin belt-bounded, radar-smooth fill),
coronae and Lakshmi Planum excluded, and incompleteness below ~10^5 km².
5. **Byrne 2021 (PNAS) uses the same machinery as `planet_tectonics.mantle_flow`**: radial-viscosity flow kernels
   driven by density loads, giving lithosphere stresses. On a Venus-type planet the same code path, with a thin
   brittle layer over a weak lower crust (5–20 km instead of Earth's ~50 km), is how blocks would break up and move
   without subduction. This is the strongest hint that the existing model could reach a Venus regime by changing
   inputs rather than adding a new system.
