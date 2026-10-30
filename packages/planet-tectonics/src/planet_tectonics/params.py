"""The plate layer's parameters (PARAMS: every tunable, with its source and history), the boundary-type codes and
Earth's stable-crust thickness quantiles."""


PARAMS = {
    'dt_myr': 1.,
    'ocean_thick_m': 7000.,
    'cont_ref_thick_m': 38000.,    # continental crust standing at cont_ref_elev_m
    'cont_ref_elev_m': 300.,
    'rho_crust': 2830., 'rho_mantle': 3300., 'rho_water': 1030., 'rho_sediment': 2400.,
    # sediment compaction (item B-a, thirty-sixth pass). With it, `sed` is the solid (grain) thickness, of density
    # rho_grain, and the elevation reads the bulk thickness and mass of a column whose porosity falls with depth below
    # its top by Athy's law, phi0 exp(-c z) (Sclater & Christie 1980: shale phi0 0.63, c 0.51 /km on oceanic crust;
    # shaly sandstone 0.56, 0.39 /km on continental crust); pores are water-filled (below sea level, and as ground
    # water in a dry column). A 0.5 km deep-sea cover (bulk) then stands at ~1,750 kg/m^3, a 5 km basin at ~2,250-2,300. Every
    # ledger stays in solid volume (eroded crust becomes grain volume 1:1; the grains are 6% lighter than the crust).
    # Within a plate step the erosion kernels move solid heights (graph_landscape unchanged, so the 11.53 km stage is
    # untouched); the next elevation call reads bulk. False: `sed` is a layer of rho_sediment (all runs before 26 Oct).
    # Default True since the thirty-seventh pass (sweep-20261025-201623, 3 seeds: gates pass, guards hold; deep floor
    # 4,317 m vs Earth 4,205 +- 150, the control 4,524; ocean bulk sediment 569 m vs GlobSed 610)
    'sediment_compaction': True, 'rho_grain': 2650.,
    'phi0_ocean': .63, 'compaction_c_ocean_per_km': .51, 'phi0_continent': .56, 'compaction_c_continent_per_km': .39,
    # sea level from a fixed water inventory (Earth's ocean, Charette & Smith 2010): each elevation is returned relative
    # to the level that holds this volume, so a deeper or wider ocean lowers the sea instead of drowning continents.
    # Elevations are first computed on the fixed reference datum (ridges at -2,600 m, cont_ref_elev_m). None = sea
    # level fixed at the datum (before step 8)
    'water_volume_km3': 1.335e9,
    'te_m': 30e3,                  # elastic thickness for the elevation (flexural isostasy of crust and sediment)
    # continental crust thinner than cont_ref_thick_m stands at its own (local, Airy) depth instead of loading the
    # elastic plate with a negative load (eighteenth pass). Crust thinned by stretching is compensated locally (the
    # McKenzie 1978 stretching model assumes local isostasy; young rifted lithosphere is weak, Te of a few km, Watts
    # 2001). As a flexural load a narrow strip thinned to 10-20 km lifted the young ocean beside it to +2-12 km (Te 30
    # km spreads its rebound over ~70 km) and sat 2.5 km below its own depth. Thickening and the erosion of thick belts
    # stay flexural. On since 24 Oct (sweep-20261024-114002: passes the elevation tests on 3 seeds, 0% of young
    # coastal ocean above sea level vs 31-33% without)
    'thinning_local_isostasy': True,
    # ocean floor eroded below its sediment loses the rest from its crust (eighteenth pass): before, the clamp
    # dsed >= -sed dropped it from the ledger while the material was deposited elsewhere (8 km^3/yr created from
    # nothing with Earth rifting). False = the old leak, for reproducing earlier runs only
    'ocean_crust_erosion': True,
    # sediment on continental crust as its own layer (twentieth pass): deposits on continental nodes go to `sed`
    # (rho_sediment, a flexural load like ocean-floor sediment: dry where its top stands above sea level, displacing
    # water below), and erosion takes that sediment first, then the crust. Before, deposits became crust of rho_crust:
    # filling 1 m of water on a submerged margin took (rm-rw)/(rm-rc) = 4.8 m of material instead of
    # (rm-rw)/(rm-rs) = 2.5 m, so margins trapped ~2x Earth's volume and old rifted margins thickened into abrupt
    # edges (nineteenth pass, margin_age_probe). Sediment on continental crust consumed in a collision joins the
    # receiving belt's crust (sediments are folded and metamorphosed into orogens). CRUST1.0's total thickness counts
    # crust + sediment, so the reviews (crust_review, margin_width) compare that sum. False = deposits booked as crust
    # (all runs before 24 Oct). On since 24 Oct (sweep-20261024-143735: passes the elevation tests on 3 seeds)
    'continental_sediment': True,
    'threads': None,              # threads for the spherical-harmonic flexure (None: all cores)
    # coarse erosion: the landscape model's surface processes at the plate grid (calibrated: scripts/calibrate_coarse_erosion)
    'erosion': True, 'erosion_substeps': 4,
    # K factor, incision threshold theta (A^m S at 46 km) and slope cap Sc from Earth: the law fitted to OCTOPUS 10Be
    # basin-averaged erosion rates on ETOPO at this grid (fifteenth pass, scripts/erosion_law_check.py). Every sweep
    # passed these explicitly from 23 Oct; defaults since 29 Oct. Before: K 0.75 (sweep-20261018-164500), theta 8
    # (step 6) then 6 (sweep-20261022-194148), Sc 0.05; reviews read old runs with crust_review.LEGACY_DEFAULTS
    'erosion_K_factor': 1.0, 'erosion_Sc': .036, 'erosion_Kd_factor': 4.,
    'erosion_theta': 200.,
    'stable_erodibility': .15,     # old crust relative to fresh orogens (as the test bed)
    'orogen_erodibility_tau_myr': 50.,
    # which crust resets the orogeny clock (erodibility 1, decaying with orogen_erodibility_tau_myr) (twenty-fourth
    # pass). The 46 km erosion law was calibrated with erodibility 1 inside Bird's (2003) orogens (16% of Earth's land,
    # 76% of its land above 2 km; scripts/erosion_law_check). 'collision' (until now): only nodes receiving collision
    # volume, which left shortening belts and collapsed crust at the stable erodibility (young-orogen erodibility on
    # 1-2% of the planet's land, twenty-third pass). 'belt': also the subduction belt on overriding continental crust
    # (within 2 sigma of subduction_d0_cont_km, where 95% of the belt's added volume goes, while the trench
    # converges) and crust in the collapse regime (thicker than collapse_thick_m: hot, weak crust). 'strain': also
    # every continental node whose crust thickens or thins tectonically faster than orogen_strain_rate
    # (|tectonic thickening| / thickness: the areal dilatation of incompressible crust). 3.5e-9 /yr is measured:
    # the threshold at which GSRM v2.1's areal dilatation (Kreemer et al. 2014) flags as much of Earth's land as
    # Bird's orogens cover (17%; it flags 59% of orogen land, 8% of the rest, 67% of land above 2 km;
    # scripts/earth_deformation). Default 'belt' since the interior sweep (sweep-20261024-193840, twenty-fifth pass;
    # 'collision' before)
    'orogen_rule': 'belt', 'orogen_strain_rate': 3.5e-9,
    # subduction zones (overriding plate): thickening flux per unit trench length = convergence rate x H, spread as a
    # Gaussian band at distance d0 from the trench; oceanic overriding plates: island arcs, which become continental
    # crust once thick enough.
    # Continental margins, two modes:
    #  * margin_shortening False (the 18 Oct Andes-type belts, sweep-20261018-144253): H_cont 2,500 m of *new* crust per
    #    m of convergence. It stood in for Andean shortening + arc magmatism but added 5-8 km^3/yr of volume (Earth's
    #    arc addition at continental margins ~1.0; Scholl & von Huene 2009), which the trenches recycled (step 8,
    #    follow-up 1, tenth pass).
    #  * margin_shortening True (23 Oct): upper-plate shortening moves existing crust into the belt, volume conserved
    #    per trench segment: rate s = max(0, upper plate's trenchward speed - rollback_mm_yr) (Oncken et al. 2006:
    #    Central Andes shortening ~ upper-plate westward drift - slab rollback, which held at ~17 mm/yr; present
    #    shortening 9-12 mm/yr). The volume s x T (T: local crust thickness) per unit trench length comes from a
    #    foreland band at foreland_d0 +- foreland_sigma behind the belt (underthrust foreland; lateral crust flow
    #    refills it from the interior; the plate raster cannot lose area inside a rigid plate) and is added to the
    #    belt at d0 +- sigma. New volume is only arc magmatism, arc_H_cont_m x convergence (750 m: ~40 km^3 per km of
    #    trench per Myr at 54 mm/yr, Earth's ~1.0 km^3/yr over ~25,000 km of continental-margin trench)
    'subduction_H_cont_m': 2500., 'subduction_d0_cont_km': 150., 'subduction_sigma_cont_km': 60.,
    # Every sweep passed True from 23 Oct (Earth's margin budget, thirteenth and sixteenth passes); default since 29 Oct
    # (before: False, with subduction_H_cont_m as new crust)
    'margin_shortening': True, 'rollback_mm_yr': 17., 'arc_H_cont_m': 750.,
    'foreland_d0_km': 350., 'foreland_sigma_km': 100.,
    'arc_H_ocean_m': 300., 'arc_d0_km': 150., 'arc_sigma_km': 50.,
    'arc_to_continent_thick_m': 20000.,
    # B-d (thirty-ninth pass): arc crust docks instead of subducting. With arc_accretion, oceanic crust thicker than
    # ocean_thick_m that a trench consumes hands its excess to the overriding plate's node as volume (as collision
    # volume): all of it where the column is at least arc_buoyant_thick_m (Cloos 1993: island arcs thicker than ~15 km
    # are buoyant enough to collide), else its upper arc_accrete_upper_m (Tetreault & Buiter 2014, Solid Earth 5: "only
    # the upper 2-5 km of arc crust are accreted onto continents through thin-skinned thrusting"); the normal 7 km of
    # basement subducts. Received on continental crust it is an arc-continent collision (orogeny marked). Default off
    # until a sweep decides.
    'arc_accretion': False, 'arc_buoyant_thick_m': 15000., 'arc_accrete_upper_m': 3500.,
    # B-d2 (fortieth pass): what feeds the trench fluxes (continental and island arcs, subduction erosion). 'band': each
    # overriding trench node adds closing rate x H x a Gaussian in trench distance, so the flux follows the node-typed
    # trench, which on the planet's rough boundaries runs 1.2-2.7x the true line (ocean arcs 1.7-2.7x the contour's
    # intra-oceanic consumption; on three of six seeds more arc crust than the whole planet's subduction consumed).
    # 'consumed': the oceanic area a trench consumed last step (m^2, at the world node it disappeared under) goes to the
    # nearest overriding trench node of that plate; each term adds area x H, spread over that node's band with its
    # Gaussian, normalized over the band (so a straight trench gets the 'band' law with v = consumed area / length).
    # Earth's H values are per m of convergence per m of trench, i.e. per m^2 consumed. Default 'band' until a sweep decides.
    'trench_flux': 'band',
    # B-d2 (fortieth pass): which side of an ocean-ocean subduction boundary is the trench (downgoing) side. 'node': each
    # node by its own crust age against its neighbour's (as built); the overlap resolution meanwhile consumes crust by one
    # polarity per plate pair (_ocean_polarity), and the two disagree at ~22% of oceanic overriding-typed nodes (seed 1001,
    # -245..-178 Myr): arcs grow on the downgoing plate and slab pull goes to the overriding one there. 'pair': where the
    # pair overlapped this step, the node follows that polarity; elsewhere the node rule. Default 'node' until a sweep decides.
    'trench_polarity': 'node',
    # subduction erosion: the trench removes crust from the base of the overriding plate near it (forearc), the
    # main continental loss balancing arc growth on Earth (~1.3-1.7 km^3/yr); flux = convergence x H_se. On oceanic
    # crust only arc crust above ocean_thick_m can be removed
    'subduction_erosion_H_m': 700., 'subduction_erosion_d0_km': 50., 'subduction_erosion_sigma_km': 40.,
    # share of trench sediment scraped onto the overriding plate. Earth (von Huene & Scholl 1991): ~20% at small and
    # medium prisms (16,300 km of trench), ~30% at large prisms (8,200 km), none at nonaccretionary margins
    # (19,000 km): 0.13 length-weighted. 0.75 (18 Oct) matched Earth's *absolute* sediment subduction (~1.65
    # km^3/yr) while the deep-sea sink deleted half the export; with the sink closed it returned three quarters of
    # all exported sediment to the continents (a loop of gain up to 4; tenth pass). Every sweep passed 0.13 from 23 Oct;
    # default since 29 Oct (before: 0.75)
    'sediment_accretion': .13,
    # continental rifting: a gap opening between continental crust is stretched continental crust (rift_thinning x its
    # neighbours' thickness, volume taken from them) until it would be thinner than breakup_thick_m; then ocean floor
    # breakup 17.5 km since 22 Oct (sweep-20261022-194148; 20 km left rifted margins one cell wide and only ~15% of the
    # continental crust under water, Earth ~29%)
    # 10 km since 24 Oct (twenty-first pass, with rift_rule pure_shear; Earth's hyper-extended margins break up at
    # ~5-15 km); 17.5 km was fitted with the node rule (sweep-20261022-194148)
    'rift_thinning': .6, 'breakup_thick_m': 10000.,
    # the volume for a rift gap comes from continental crust within rift_width_km of it (linear taper), not only
    # from its direct neighbours, so rifted margins thin over a necking zone (Earth: ~50-200 km) instead of one cell.
    # 100 km since 24 Oct (the realistic-margin value used in every sweep since 23 Oct; 0 before)
    'rift_width_km': 100.,
    # how a rift gap's crust is set (twenty-first pass). 'node': rift_thinning x its continental neighbours' mean
    # thickness, a factor per grid cell (38 -> 22.8 -> 13.7 km -> breakup: ~2 cells of stretched crust per side at any
    # resolution, 45-70 km of thin crust per km of passive edge against Earth's 140). 'pure_shear': the gap is the new
    # area of a stretching zone of half-width rift_width_km (linear taper of the strain rate from the gap cell's edge,
    # through continental crust only); every zone node thins as pure shear, T exp(-c) with c its share of the opened
    # area, and the gap receives exactly the volume the zone loses, i.e. the strain-weighted mean thickness of the
    # zone. The gap breaks up (ocean floor) when that mean is below breakup_thick_m. Stated in km, so a margin's profile
    # does not depend on the cell size; a stretched zone keeps thinning while the rift stays open, which gives the
    # hyper-extended domain (crust ~10 km over 70-200 km, Brune et al. 2014) and a necking zone at its inner edge.
    # pure_shear since 24 Oct (sweep-20261024-173200, 3 seeds: passes the elevation tests; thin crust per km of passive
    # edge 130-144 km at 1 degree vs 58-83 with 'node', Earth 140). 'node' reproduces earlier runs
    'rift_rule': 'pure_shear',
    # pelagic (biogenic and wind-blown) sediment on the ocean floor, m/Myr (Earth's deep sea ~1-10; GlobSed mean
    # thickness ~0.5 km); without it the planet's ocean floor carried ~0.1 km of sediment and stood 150-300 m too deep
    'pelagic_m_myr': 0.,
    # continental shelves (planet_tectonics.shelf): wave transport in shallow water, kappa_wave * exp(-depth / wave_base),
    # and coastal retreat at cliff_retreat_m_yr damped by the water in front of the coast; run after every coarse
    # erosion substep. shelf_kappa_m2_yr = 0 turns it off. On since 22 Oct (sweep-20261022-194148: 3e4 / 0.003 with
    # breakup 17.5 km; on one-cell margins it had no effect, sweep-20261021-182714)
    'shelf_kappa_m2_yr': 3e4, 'shelf_wave_base_m': 50., 'cliff_retreat_m_yr': .003,
    # deep-sea sediment: the landscape kernels' marine step places sediment only on shelves (within 150 km of land) and
    # drops the rest ("to the deep sea", right for a regional window). On the planet that deleted ~7 km^3/yr, more
    # than the continents' net loss (step 8 follow-up 1). Kept, the rest ponds where its path ended and fills the
    # sea-floor lows from the bottom up (graph_landscape.fill_and_spill: turbidites, abyssal plains), then rides the
    # ocean floor to the trenches like any other sediment. (A first version spread it by diffusion over the connected
    # open sea; enclosed pockets could not pass it on and grew 20-250 km towers.) False = the old sink
    'keep_deep_sediment': True,
    'mantle': True,                # dynamic topography and hotspots (planet_tectonics.mantle)
    # dynamic topography (planet_tectonics.mantle): std 800 m to spherical-harmonic degree 30 with a flat (slope 1) spectrum,
    # chosen by the step 8 sweep (sweep-20261017-223508: texture 0.23 vs 0.77 for the old 400 m / degree 12 / slope 2.5;
    # residual topography on Earth is +-~1 km with power to degree ~30, Hoggard et al. 2016, the high end of estimates)
    # (sweep-20261018-144253, scored with arrangement and budget too, then preferred 500 m)
    'dyn_std_m': 500., 'dyn_lmax': 30, 'dyn_slope': 1., 'dt_knot_myr': 40.,
    # plate splitting (Poisson events). Rifting: rate per Myr for a 100 Mkm^2 plate, scaled by area and doubled for
    # mostly continental plates; the split line prefers old sutures; halves diverge at split_speed_cm_yr.
    'split_rate_per_myr': .006, 'split_min_area_mkm2': 12., 'split_speed_cm_yr': 2., 'suture_weakness': .6,
    # subduction initiation: ocean older than init_age_myr on a plate breaks away and converges on the rest
    'init_age_myr': 170., 'init_rate_per_myr': .004, 'init_speed_cm_yr': 3.,
    'init_area_ref_mkm2': 2.,      # initiation rate scales with old-ocean area / this (capped at 5x)
    # suturing: two plates in continental collision along more than merge_contact_km for merge_after_myr in a row
    # become one plate (as India is joining Eurasia); the smaller plate's crust moves into the larger plate's frame
    'merge_contact_km': 1000., 'merge_after_myr': 40.,
    'merge_max_ratio': .35, 'merge_max_area_mkm2': 150.,   # docking: a smaller block onto a bigger one, size-capped
    # plate pieces (item B-b, thirty-seventh pass). A plate's world territory can fall into separate pieces (a plate
    # growing inside another cuts its host in two; overrides at trenches); with 'shared' (as built) every piece keeps
    # the plate's one Euler pole, thousands of km away: 6-13% of the surface at 0 Myr, 35-53% of all boundary and half
    # the intra-oceanic trench (thirty-sixth pass). 'plates': on Earth a plate is one connected rigid piece of
    # lithosphere. A detached piece of at least piece_min_area_mkm2 becomes its own plate (the parent's Euler vector at
    # the split, then its own forces); a smaller one joins the plate around most of its edge, which keeps its motion
    # (microplate capture: the Monterey, Arguello and Magdalena microplates by the Pacific plate, Nicholson et al.
    # 1994); a piece with no interior node (a staircase crumb on a boundary, 1-5 nodes, 0.07% of the surface) is left.
    # 1 Mkm^2: Bird's (2003) 20 plates of at least 1 Mkm^2 hold 98% of Earth's surface; the 32 smaller ones (10 Mkm^2)
    # live in deforming zones the n = 217 grid barely resolves. Default 'plates' since the thirty-eighth pass
    # (sweep-20261025-220434: detached max 0.0018-0.0033 of the surface vs 0.06-0.095; volume -10.8% vs -12.5%)
    'piece_rule': 'plates', 'piece_min_area_mkm2': 1.,
    # sediment-filled ocean floor becomes continental crust (margin prisms, filled remnant basins), keeping its
    # elevation: oceanic nodes with more than basin_to_continent_sed_m of sediment, or any sediment above sea level
    'basin_to_continent_sed_m': 10000.,
    # off since step 8: the conversion set the new continental thickness to match the node's elevation, not its mass
    # (1 km of sediment could become ~30 km of crust); over 250 Myr it added continent over 12.8% of the surface,
    # thinned the mean continental crust from 38.5 to 34.9 km and hid a ~30% continental volume loss. Thick fills now
    # stay sediment on oceanic crust, as in the South Caspian or the Gulf of Mexico
    'basin_to_continent': False,
    # lateral spreading of continental crust (m^2/yr): kappa_normal everywhere, plus gravitational collapse of thick
    # (hot, weak) crust, crust_kappa * (exp((T - collapse_thick) / collapse_scale) - 1) up to collapse_max_factor x
    # crust_kappa, which keeps plateaus near ~70 km as Tibet and the Altiplano. Cold crust of normal thickness barely
    # flows: kappa_normal = 1000 erased all 100 km crustal structure within ~10 Myr (step 8 texture check).
    # kappa_normal 30 since 24 Oct (twentieth pass, sweep-20261024-143735: passes the elevation tests on 3 seeds):
    # kappa ~ rho_c g (1 - rho_c/rho_m) h^3 / (3 eta) for a 10 km lower-crust channel is ~40 m^2/yr at eta = 1e21 Pa s
    # (stable continental lower crust; 1000 is ~4e19, hot crust), and Earth's 180 Myr Atlantic margins survive,
    # which needs kappa << 125. At 1000 the flow thickened old passive margins by 5-10 km per 100 Myr
    'crust_kappa_normal_m2_yr': 30.,
    'crust_kappa_m2_yr': 1000.,
    'collapse_thick_m': 45000.,
    'collapse_scale_m': 6000.,
    'collapse_max_factor': 200.,
    # collapse past the continent edge (seventeenth pass): the collapse part of the crust flow (crust thicker than
    # collapse_thick_m) may also cross the coastline into oceanic nodes of the same plate that are not on a plate
    # boundary, volume-conserving; such nodes become continental crust through the arc rule (arc_to_continent_thick_m).
    # On Earth thick orogens collapse and extend, adding continental area (Basin and Range, Aegean); before this, crust
    # flow stopped at the coastline, so collisions only ever removed continental area. Normal-thickness crust does not
    # cross (cold crust is strong). Default False until the 3-seed sweep decides
    'collapse_edge': False,
    # thickness of the supercontinent's crust at the start (twenty-fourth pass). 'assumed' (until now): 37.5 km + 2 km
    # beyond 300 km from old sutures + 4 km near them + a +-1.5 km field (interior p50 39.2-40.1 km, p10-p90 3-4 km
    # wide); it survives 250 Myr almost untouched (95-99% of the present interior, twenty-third pass). 'earth': the
    # same field quantile-mapped onto Earth's stable continental crust (EARTH_STABLE_CRUST_KM: CRUST1.0 crystalline
    # thickness of Precambrian and platform crust types more than 500 km from plate boundaries, p10 / 50 / 90 =
    # 32.9 / 38.2 / 43.0 km), keeping its spatial structure (thicker near old sutures); the thinned rim at the coast
    # as before. Default 'earth' since the interior sweep (sweep-20261024-193840, twenty-fifth pass; 'assumed' before)
    'start_crust': 'earth',
    # the start's plate geometry (thirty-eighth pass). 'voronoi' (as built): plates are Voronoi cells of the whole
    # sphere, n_cont_plates seeded inside the supercontinent, so its plates reach into the ocean and its coast is a
    # passive margin inside plates (continental-side trench at -250 Myr 0.2-0.7k km); a closed sphere consumes what it
    # opens, so all closing is ocean-ocean (intra-oceanic trench 33-38k km at the start, 0.40 of the ocean-ocean
    # boundary; Earth 20.8k, 0.14; thirty-seventh pass). 'coast': the same seeds grow only within their own domain
    # (continental plates end at the coast, oceanic plates fill the ocean), so the supercontinent's rim is a plate
    # boundary all round and its type follows from the motions (the continental plates move away from the centre, as
    # before): converging stretches subduct ocean under the continent, diverging ones open new ocean that joins the
    # continental plate (a passive margin). Pangea was ringed by subduction along its Panthalassa side (the Terra
    # Australis orogen, ~18,000 km, Cawood 2005; the circum-Pangea subduction girdle, Collins 2003), with passive
    # margins mainly on the opening Tethys side. An ocean (or continent) component without a seed gets its own plate.
    # No new number
    'start_plates': 'voronoi',
    # the oceanic start plates' motion (thirty-eighth pass). 'random' (as built): random poles, 3-7 cm/yr; the force
    # law turns them over tau_axis_myr. 'forces': each oceanic start plate gets the speed and pole the plate layer's own
    # force law (_driving, fully relaxed) gives on the start's geometry, with the slabs Pangea had: ocean nodes on the
    # coast are trench (the subduction ring; meaningful with start_plates 'coast') and ocean-ocean boundaries ridges
    # (as the start's age field has them). Continental plates keep their dispersal motion. No new number
    'start_poles': 'random',
    'v_min_cm_yr': 1.5, 'v_slab_cm_yr': 9.,  # plate speed without / with a full slab
    'slab_ref_km': 6000.,          # trench length giving the full slab-pull speed
    'v_ridge_cm_yr': 1.5,          # extra speed from ridge push with ridge_ref_km of ridge
    'ridge_ref_km': 8000.,
    'ridge_torque': .3,            # weight of a ridge node relative to a trench node in turning the pole
    'cont_drag': .3,               # speed factor (1 - cont_drag * continental fraction)
    'collision_ref_km': 1500.,     # collision length that halves the target speed
    'tau_speed_myr': 10., 'tau_axis_myr': 20.,
    'boundary_min_mm_yr': 5.,      # below this, a boundary is a transform
    # boundary normal used to type boundaries (thirty-third pass). 'node' (as built): each node is typed by the one
    # neighbour on another plate with the largest normal rate, so the normal is a grid direction; on the grid's
    # staircase a transform that does not follow the grid shows normal rates of ~0.7 v against ~0.7 v of shear and is
    # typed ridge on one step and collision or subduction on the next (pure transforms on 12 random great circles:
    # median 30% of nodes typed transform, min 0%, at n = 48 and 217; a run's transforms 11% of node-typed boundary vs
    # 36% of contour-typed length, scripts/boundary_lengths). 'ring': the normal toward a neighbouring plate b is the
    # tangent direction of the sum over the node's two-ring of (x_k - x_i) (+1 on b, -1 on its own plate), the local
    # line's normal (96-100% typed transform on the same test). No new number. Default 'ring' since the thirty-fourth
    # pass (sweep-20261025-170555: its control passes the gates and guards on 3 seeds, node typing failed the elevation
    # tests on 2 of 3 with young-ocean offsets -356 / -307 m); 'node' reproduces the runs before it step by step
    'boundary_normal': 'ring',
    'hidden_steps': 3,             # overlapped crust is consumed only after this many steps under another plate
}



RIDGE, SUBDUCTION, COLLISION, TRANSFORM = 1, 2, 3, 4



# Earth's stable continental crust (start_crust 'earth'): crystalline thickness (CRUST1.0 crsthk - sedthk, Laske et
# al. 2013) of Precambrian and platform crust types (D-, E-, F-, G1, G2, H1, H2, I1, I2) more than 500 km from Bird's
# (2003) plate boundaries, area-weighted on the n = 217 graph, at the quantiles EARTH_STABLE_CRUST_Q (twenty-fourth
# pass; tests/test_plates.py recomputes it)
EARTH_STABLE_CRUST_Q = (.01, .05, .10, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90, .95, .99)


EARTH_STABLE_CRUST_Q = (.01, .05, .10, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60, .65, .70, .75, .80, .85, .90, .95, .99)
EARTH_STABLE_CRUST_KM = (25.65, 30.54, 32.91, 34.36, 35.2, 35.9, 36.43, 36.9, 37.19, 37.7, 38.17, 38.7, 38.99, 39.44,
                         39.99, 40.51, 41.02, 42.0, 43.0, 44.9, 48.79)
