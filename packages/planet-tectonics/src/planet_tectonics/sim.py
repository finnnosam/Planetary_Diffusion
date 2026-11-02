"""Plate layer, part 1 (docs/planet-layer-plan.md, step 5): rigid plates on the cube sphere, each keeping its crust in
its own rotating frame; ridges, subduction, collision thickening and ocean-floor ages; a supercontinent start.

Representation
--------------
Every plate owns a full raster on the same CubeSphere node set, interpreted in the plate's own frame, plus a
rotation Q (plate frame -> world). Each step (default 1 Myr):

1. rotate every plate about its Euler pole;
2. world node p is *claimed* by plate k if k's raster is set at the node nearest Q_k^T p;
3. no claim = divergence: new oceanic crust (7 km, age 0) joins the nearest plate (ridge);
   several claims = convergence: continental crust beats oceanic, younger ocean beats older; the losing ocean is
   subducted; losing continental crust under continental crust is added to the winner's thickness (collision);
4. thick crust spreads laterally (diffusion of thickness within continents: gravitational spreading of belts);
5. each plate's raster is updated by *forward* lookup (plate node j -> world node nearest Q_k j): nodes the plate
   still owns keep their values and receive only *increments* (collision, spreading); newly owned nodes are
   initialised from the world; nodes it lost are cleared. Crust is never resampled from itself, so 250 Myr of
   motion adds no interpolation blur.
6. plate motion: slab pull turns each plate's pole so that it moves toward its trenches and speeds it up with
   trench length (Forsyth & Uyeda 1975); continental area and collisions slow it.

Elevation here is isostatic only (ocean depth from age, GDH1; continents from Airy crustal thickness with water
loading below 0 m). Erosion, rifting, terrane transfer, hotspots and mantle swells are part 2 (step 6).
Times are in Myr (negative = before present), distances in m, speeds in cm/yr.

The simulator class; its methods are grouped by mechanism in start, composition, crust, surface, motion and events."""
import numpy as np

from planet_grid import fields
from planet_grid.sphere_flexure import SphereFlexure
from planet_landscape import graph_landscape as gl

from planet_tectonics.crust_kernels import edge_spread
from planet_tectonics.isostasy import sediment_column
from planet_tectonics.lithosphere import stretched_age
from planet_tectonics.mantle import Mantle
from planet_tectonics.params import PARAMS
from planet_tectonics.start import StartMixin
from planet_tectonics.composition import CompositionMixin
from planet_tectonics.crust import CrustMixin
from planet_tectonics.surface import SurfaceMixin
from planet_tectonics.motion import MotionMixin
from planet_tectonics.events import EventsMixin


class PlateSim(StartMixin, CompositionMixin, CrustMixin, SurfaceMixin, MotionMixin, EventsMixin):
    """Plates on a planet_grid.sphere_graph.CubeSphere (n = 217 for 46 km)."""


    def __init__(self, graph, seed=1, n_plates=12, n_cont_plates=5, cont_fraction=.38, t0=-250., start='supercontinent',
                 params=None):
        self.g = graph
        self.p = {**PARAMS, **(params or {})}
        self.t = float(t0)
        self.rng = np.random.default_rng(seed)
        self.bound = graph.diffusion_bound()
        self.h_min = graph.spacing_min()
        # loads are passed as mass per area (kg/m^2): subaerial loads displace mantle against air (rho_m), submarine
        # loads replace water and the depression fills with water (rho_m - rho_w)
        self.flex_dry = SphereFlexure.for_graph(graph, self.p['te_m'], rho_load=1., rho_mantle=self.p['rho_mantle'],
                                                nthreads=self.p['threads'])
        self.flex_wet = SphereFlexure.for_graph(graph, self.p['te_m'], rho_load=1.,
                                                rho_mantle=self.p['rho_mantle']-self.p['rho_water'], nthreads=self.p['threads'])
        mp = {k: self.p[k] for k in ('dyn_std_m', 'dyn_lmax', 'dyn_slope', 'dt_knot_myr')}
        if self.p['dyn_rule'] != 'random':           # D: the slab flow (planet_tectonics.mantle_flow)
            mp.update({k: self.p[k] for k in ('dyn_rule', 'rho_mantle', 'rho_water', 'threads', 'slab_dip_deg',
                                              'lm_sink_mm_yr', 'flow_z_cut_km', 'flow_lmax', 'viscosity_file',
                                              'slab_alpha')})
        self.mantle = Mantle(graph, np.random.default_rng([seed, 7]), t0, params=mp) if self.p['mantle'] else None
        self.dyn = self.mantle.topography(t0) if self.mantle else np.zeros(graph.size)
        self.next_id = 0
        self.plates = []
        self.events = []
        self.ledger = dict.fromkeys(('new_nodes', 'holes', 'increments', 'collision_out', 'cont_subducted', 'cont_lost_other'), 0.)
        # sediment on continental crust (continental_sediment_volume), by cause; with continental_sediment its
        # change over a run equals the sum of these (tests/test_plates.py)
        self.sed_ledger = dict.fromkeys(('increments', 'arc_continent', 'new_nodes', 'holes', 'collision_out',
                                         'cont_subducted', 'merge'), 0.)
        self.erosion_totals = {'continental_eroded_m3': 0., 'continental_deposited_m3': 0., 'ocean_sediment_m3': 0.}
        # continental area (flag) gained and lost by mechanism, m^2 (the area ledger; the rest of the net change is
        # nearest-node jitter, holes and merges, read as the residual)
        self.area_ledger = dict.fromkeys(('collision', 'subduction', 'rift', 'arc'), 0.)
        if start != 'supercontinent':
            raise NotImplementedError('only the supercontinent start exists so far (dispersed start: later)')
        self._init_supercontinent(n_plates, n_cont_plates, cont_fraction)
        self.compose_from_plates()
        if self.p['start_poles'] == 'forces':
            self._start_poles_from_forces()
        elif self.p['start_poles'] != 'random':
            raise ValueError(f"start_poles {self.p['start_poles']!r}: 'random' or 'forces'")
        self.area0 = float(graph.area[self.world['cont']].sum())


    # ------------------------------------------------------------------ stepping
    def step(self):
        g, p = self.g, self.p
        dt = p['dt_myr']
        for pl in self.plates:
            pl.Q = fields.rotation(pl.axis, pl.speed*1e-2*dt*1e6/g.radius)@pl.Q
        self.t += dt
        if self.mantle is not None:
            self.mantle.advance(self.t, dt)
            if self.mantle.flow is not None:
                self._feed_slab_record(dt)
                wet = self.z_last < 0 if hasattr(self, 'z_last') else None
                self.dyn = self.mantle.topography(self.t, wet)
            else:
                self.dyn = self.mantle.topography(self.t)
        idx, claimed = self._claims()
        self._world_from(idx, claimed, record=True)
        w = self.world
        rift_dthick = self._rift()
        # lateral spreading of thick crust inside continents (conservative; increments only)
        before = w['thick'].copy()
        th = before
        T = dt*1e6
        while T > 0:                                     # nonlinear: substeps follow the current largest diffusivity
            kap = p['crust_kappa_normal_m2_yr']+p['crust_kappa_m2_yr']*np.minimum(
                np.expm1(np.clip(th-p['collapse_thick_m'], 0, None)/p['collapse_scale_m']), p['collapse_max_factor'])
            kmax = kap[w['cont']].max() if w['cont'].any() else 0.
            if kmax <= 0:
                break
            h = min(T, .9/(kmax*self.bound))
            th = gl.diffuse_variable(th, w['cont'], kap, h, g.nbr, g.dist, g.axis, g.area)
            if p['collapse_edge'] and hasattr(self, 'btype'):
                kc = p['crust_kappa_m2_yr']*np.minimum(np.expm1(np.clip(th-p['collapse_thick_m'], 0, None)/p['collapse_scale_m']),
                                                     p['collapse_max_factor'])
                th1 = edge_spread(th, w['cont'], kc, np.clip(th-p['collapse_thick_m'], 0, None), self.owner,
                                  (~w['cont']) & (self.btype == 0), h, g.nbr, g.dist, g.axis, g.area)
                self.erosion_totals['edge_collapse_m3'] = self.erosion_totals.get('edge_collapse_m3', 0.)+float(
                    (np.maximum(th1-th, 0)*g.area)[~w['cont']].sum())
                th = th1
            T -= h
        dthick = th-before+rift_dthick
        thin = w['cont'] & (rift_dthick < 0) & (before > 0)
        if thin.any():                                   # C: uniform stretching thins the lithosphere with the crust
            ratio = np.clip((before[thin]+rift_dthick[thin])/before[thin], 0., 1.)
            w['therm'][thin] = self.t-stretched_age(self.t-w['therm'][thin], ratio, p)
            self.therm_changed |= thin
        th = th+rift_dthick
        dsed = np.zeros(g.size)
        sub = self._subduction_thickening()
        if self.mantle is not None:
            vol = self.mantle.volcanic(self.t, dt)
            self.erosion_totals['volcanic_added_m3'] = self.erosion_totals.get('volcanic_added_m3', 0.)+float((vol*g.area).sum())
            sub = sub+vol
        dthick = dthick+sub
        th = th+sub
        self.tectonic = dthick.copy()                    # tectonic thickening this step, m (no erosion); + collisions below
        arc = (~w['cont']) & (th >= p['arc_to_continent_thick_m'])
        if arc.any():                                    # island arcs grown thick enough are continental crust
            self.area_ledger['arc'] += float(g.area[arc].sum())
            w['cont'] = w['cont'] | arc
            w['therm'][arc] = self.t                     # C: new continental crust, new lithosphere
            self.therm_changed |= arc
            self.new_continent = arc
        else:
            self.new_continent = np.zeros(g.size, bool)
        if p['erosion']:
            e_cont, dsed = self._erode(th)
            dthick = dthick+e_cont
            th = th+e_cont
        w['thick'] = th
        if p['pelagic_m_myr'] > 0:
            pel = np.where(w['cont'], 0., p['pelagic_m_myr']*dt)
            dsed = dsed+pel
            self.erosion_totals['pelagic_m3'] = self.erosion_totals.get('pelagic_m3', 0.)+float((pel*g.area).sum())
        w['sed'] = w['sed']+dsed
        conv_t, conv_s = self._sediment_basins_to_continent()
        dthick = dthick+conv_t
        dsed = dsed+conv_s
        self._update_rasters(dthick, idx, dsed)
        self.plates = [pl for pl in self.plates if pl.mask.any()]
        if p['piece_rule'] == 'plates':
            self._separate_pieces()
        elif p['piece_rule'] != 'shared':
            raise ValueError(f"piece_rule {p['piece_rule']!r}: 'shared' or 'plates'")
        self._classify()
        self._update_motion()
        self._merge_events()
        self._split_events()


    def _feed_slab_record(self, dt):
        """dyn_rule 'slabs': this step's parcels from the last classification's oceanic trench nodes (downgoing side):
        area closing rate x dt x BOUNDARY_NODE_H h, times the floor's thermal mass excess; they descend toward the
        overriding plate (the boundary normal at the node)."""
        from planet_tectonics.mantle_flow import thermal_mass_excess
        from planet_tectonics.motion import BOUNDARY_NODE_H
        g, p, w = self.g, self.p, self.world
        if hasattr(self, 'trench') and hasattr(self, 'b_normal'):
            rate = np.maximum(self.conv_rate, 0.)
            I = np.flatnonzero(self.trench & ~w['cont'] & (rate > 0) & (np.abs(self.b_normal).sum(1) > 0))
        else:
            I = np.zeros(0, np.int64); rate = np.zeros(g.size)
        age = self.t-dt-w['born'][I]                 # age at the last classification
        Ln = BOUNDARY_NODE_H*np.sqrt(g.area.mean())
        mass = thermal_mass_excess(age, p['rho_mantle'], p['rho_water'])*rate[I]*dt*1e6*Ln
        self.mantle.subduct(self.t, dt, g.xyz[I], self.b_normal[I] if I.size else np.zeros((0, 3)), mass, rate[I])


    def run(self, years_myr, callback=None, every=None):
        steps = int(round(years_myr/self.p['dt_myr']))
        for s in range(steps):
            self.step()
            if callback is not None and every and (s+1) % every == 0:
                callback(self)
        return self


    # ------------------------------------------------------------------ outputs
    def state(self):
        """World fields at the current time (node vectors)."""
        w = self.world
        ocean_age = np.where(w['cont'], np.nan, self.t-w['born'])
        elev = self.elevation()
        out = {'time_myr': self.t, 'plate_id': self.owner.astype(np.int32), 'continental': w['cont'],
                'crust_thickness_m': w['thick'].astype(np.float32), 'crust_age_myr': (self.t-w['born']).astype(np.float32),
                'ocean_age_myr': ocean_age.astype(np.float32), 'orogeny_age_myr': (self.t-w['orog']).astype(np.float32),
                'lithosphere_age_myr': np.where(w['cont'], self.t-w['therm'], np.nan).astype(np.float32),
                'boundary_type': self.btype if hasattr(self, 'btype') else np.zeros(self.g.size, np.int8),
                'convergence_mm_yr': (self.conv_rate*1e3).astype(np.float32) if hasattr(self, 'conv_rate') else np.zeros(self.g.size, np.float32),
                'elevation_m': elev.astype(np.float32), 'sediment_m': w['sed'].astype(np.float32),
                'sediment_bulk_m': sediment_column(w['sed'], w['cont'], self.p)[0].astype(np.float32),
                'dynamic_topography_m': self.dyn.astype(np.float32), 'new_crust': self.new_crust,
                'sea_level_m': np.float32(getattr(self, 'sea_level', 0.)),
                'tectonic_thickening_m': getattr(self, 'tectonic', np.zeros(self.g.size)).astype(np.float32)}
        fl = self.mantle.flow if self.mantle is not None else None
        if fl is not None:                           # D: the flow's coefficients and the record's profile
            rec = fl.record
            young = (rec.z >= fl.z_cut) & (self.t-rec.born < 100.)
            near = np.zeros(self.g.size, np.float32)
            if young.any():
                np.add.at(near, fields.nearest_nodes(self.g, rec.pos[young]), rec.mass[young].astype(np.float32))
            out.update({'dynamic_topography_air_m': self.mantle.topography(self.t).astype(np.float32),   # no water factor
                        'flow_lmax': np.int32(fl.lmax), 'flow_topography_alm_m': fl.H.astype(np.complex64),
                        'flow_velocity_alm_m_yr': fl.V.astype(np.complex64),
                        'slab_mass_young_kg': near,
                        'slab_mass_by_layer_kg': rec.mass_by_layer(fl.z_cut, len(fl.radii)),
                        'slab_record_kg': np.array([rec.deposited_kg, rec.total_mass(), rec.exited_kg])})
        return out


    def plate_table(self):
        g = self.g
        rows = []
        for pl in self.plates:
            own = self.owner == pl.id
            rows.append({'id': pl.id, 'area_km2': float(g.area[own].sum()/1e6), 'continental_km2': float(g.area[own & self.world['cont']].sum()/1e6),
                         'speed_cm_yr': pl.speed, 'axis': pl.axis.tolist(), 'Q': pl.Q.tolist()})
        return rows


    def continental_volume(self):
        """Continental crust volume held by the plates (m^3), including crust hidden under another plate."""
        a = self.g.area
        return float(sum((a*pl.thick*(pl.mask & pl.cont)).sum() for pl in self.plates))


    def _plate_cont_sed(self, pl):
        return float((self.g.area*pl.sed.astype(np.float64)*(pl.mask & pl.cont)).sum())


    def continental_sediment_volume(self):
        """Sediment on continental crust held by the plates (m^3), including crust hidden under another plate. Without
        continental_sediment this is only ocean-floor sediment carried onto crust that became continental (arcs)."""
        return float(sum(self._plate_cont_sed(pl) for pl in self.plates))


    def pinholes(self):
        """New-crust nodes more than 2.5 cells from any plate boundary: holes inside a plate (should be ~0). Ridge crust
        is always next to the boundary it opened at."""
        g = self.g
        nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
        edge = (self.owner[nb] != self.owner[:, None]).any(1)
        d = gl.distance_to(edge, g.nbr, g.dist, 3*np.sqrt(g.area.mean()))
        return int((self.new_crust & (d > 2.5*np.sqrt(g.area.mean()))).sum())
