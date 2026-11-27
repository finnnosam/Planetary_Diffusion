"""Pass 105 (journal 142; STATE section 0 item 1): where the excess shoreline and the missing inland lowland first
appear. Built in pass 105, read in pass 106 (an evidence checkpoint). The question and the reading rule are fixed in
journal 142 before any trace runs.

One plate-layer step, as `PlateSim.step` runs it, read as a chain of surfaces on the world grid (each its own
flexural-isostatic elevation with the sea re-solved for the fixed water volume; the probe restores `sea_level` after
every extra solve, so the traced trajectory is the untraced one: `--check-identity`):

  prev_end         the previous step's surface (what `state()` returned)
  composed_olddyn  after the plates' rotation and the plate-to-world composition, with the previous step's dynamic
                   field: the rigid move, the nearest-node composition and whatever the last return transfer left in
                   the plate rasters (nodes no world node read last step show up here)
  composed         the same with this step's dynamic field
  tectonics        after rifting, crust flow, subduction thickening, volcanism, foundering, arcs (before erosion)
  eroded           after erosion and deposition (`_erode`)
  transferred      after pelagic rain, basins turned continental and the return transfer (`_update_rasters`)
  end              after the back-arc strips and the events (`state()`)
side readings (not in the chain):
  returned         the plate rasters read back at this step's rotation right after the transfer (`transferred` as the
                   plates now hold it; every world node reads its plate node)
  end_local        `end` with local Airy isostasy in place of flexure (the flexure's share)
  end_prev_sea     `end` against the previous step's sea level (the sea solve's share)

Each surface: land fraction; one-cell holes and islands (eight-neighbour graph); unique land-water crossings per
Mkm^2; land with four dry rings and the low, flat four-ring interior (centre < 500 m, relief < 300 m), as in the
review's terrain probe; land below 200 m among land >= 300 km from water (graph distance to the nearest water node next
to land) and that far land's share; the drowned normal band (crystalline >= 30 km, -1,000 ... -200 m) and its entrants
and exits against the previous surface; the roughness (mean |value - mean of its neighbours| on continental interior
nodes) of the surface, the crystalline crust and the sediment. Holes and islands created and destroyed against the
previous surface, the created ones tallied by marks set at prev_end, before any selection on height: crystalline
thickness (< 30, 30-38, 38-45, >= 45 km), graph distance to the material continental edge (< 100, 100-300, 300-600,
>= 600 km), boundary history (within 150 km of a plate boundary at any of the last 10 traced steps).
The return transfer, directly: continental plate nodes that no world node read (area share), and the rms of the
increments they missed (the increment of the world node over them) against the rms of the applied ones.

Modes (from a warm state pickled by `warm`, so all modes start from the same state):
  default   the model as it is
  smoothed  the same after smoothing each plate's continental crust and sediment (conservative graph diffusion,
            100 km): which operation brings the roughness back, and how fast
  null      smoothed material, rigid rotation and composition only (no dynamic change, no processes, no transfer)
  null_raw  the same on the unsmoothed material

  python trace_shore.py warm  <out> --seed 1028 --stop -20           # run to -20 Myr, pickle (cloud, threads 1)
  python trace_shore.py trace <out> --mode default --steps 20
  python trace_shore.py summary <out>
  python trace_shore.py check-identity <out> --n 48 --steps 6       # traced = untraced (threads 1)
"""
import argparse
import copyreg
import importlib
import json
import pickle
import sys
import time
import types
from pathlib import Path

import numpy as np

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))

from planet_grid import fields                                            # noqa: E402
from planet_grid.sphere_graph import CubeSphere                           # noqa: E402
from planet_landscape import graph_landscape as gl                        # noqa: E402
from planet_tectonics.plates import PlateSim                              # noqa: E402

STAGES = ('prev_end', 'composed_olddyn', 'composed', 'tectonics', 'eroded', 'transferred', 'end')
SIDE = ('returned', 'end_local', 'end_prev_sea')
T_BINS = (0., 30e3, 38e3, 45e3, np.inf)
E_BINS = (0., 100e3, 300e3, 600e3, np.inf)
SMOOTH_KM = 100.


class LocalAiry:
    """Stands in for SphereFlexure: local compensation of a load (kg/m^2) against the given density."""

    def __init__(self, rho):
        self.rho = rho

    def deflection(self, load):
        return np.asarray(load, np.float64)/self.rho


class Geo:
    """Shoreline and lowland measures on one grid."""

    def __init__(self, g):
        self.g = g
        nbr = g.nbr.copy()
        bad = nbr < 0
        nbr[bad] = np.broadcast_to(np.arange(g.size)[:, None], nbr.shape)[bad]
        self.nbr8 = nbr
        i, k = np.where(g.axis & (np.arange(g.size)[:, None] < np.where(bad, -1, g.nbr)))
        self.ei, self.ej = i, g.nbr[i, k]
        self.A = g.area.sum()

    def dist_to(self, src, cutoff):
        return gl.distance_to(src, self.g.nbr, self.g.dist, cutoff)

    def masks(self, z):
        nb = self.nbr8
        land = z > 0
        hole = ~land & land[nb].all(1)
        island = land & (~land[nb]).all(1)
        return land, hole, island

    def roughness(self, v, interior):
        d = np.abs(v-v[self.nbr8].mean(1))
        a = self.g.area
        return float(np.average(d[interior], weights=a[interior])) if interior.any() else float('nan')

    def measures(self, z, T, sed, cont):
        g, a, nb = self.g, self.g.area, self.nbr8
        land, hole, island = self.masks(z)
        aland = a[land].sum()
        low, high = z.copy(), z.copy()
        for _ in range(4):
            low = low[nb].min(1); high = high[nb].max(1)
        dry4 = land & (low > 0)
        lowflat = dry4 & (high-low < 300.) & (z < 500.)
        coast_w = ~land & land[nb].any(1)
        d = self.dist_to(coast_w, 301e3)
        far = land & (d >= 300e3)
        band = (T >= 30e3) & (z >= -1000.) & (z < -200.)
        interior = cont & cont[nb].all(1)
        r = {'land': float(aland/self.A), 'holes': int(hole.sum()), 'islands': int(island.sum()),
             'crossings_per_mkm2': float((land[self.ei] != land[self.ej]).sum()/(self.A/1e12)),
             'dry4_of_land': float(a[dry4].sum()/aland), 'lowflat4_of_land': float(a[lowflat].sum()/aland),
             'low200_of_land': float(a[land & (z < 200.)].sum()/aland),
             'far300_of_land': float(a[far].sum()/aland),
             'low200_of_far300': float(a[far & (z < 200.)].sum()/max(a[far].sum(), 1.)),
             'band_of_globe': float(a[band].sum()/self.A),
             'rough_z': self.roughness(z, interior), 'rough_T': self.roughness(T, interior),
             'rough_sed': self.roughness(sed, interior)}
        return r, {'land': land, 'hole': hole, 'island': island, 'band': band, 'far_low': far & (z < 200.),
                   'lowflat': lowflat}


def pickle_save(obj, f):
    with open(f, 'wb') as fh:
        pickle.dump(obj, fh, protocol=5)


def pickle_load(f):
    with open(f, 'rb') as fh:
        return pickle.load(fh)


class Tracer:
    """Wraps one PlateSim's step operations (instance attributes shadow the methods; the step calls them as before)."""

    def __init__(self, sim):
        self.sim = sim
        self.geo = Geo(sim.g)
        self.surfaces = {}
        self.transfer = {}
        self.bhist = []
        self.prev = None
        self._wrap()

    # -- diagnostic elevation without touching the trajectory
    def zdiag(self, thick=None, sed=None, world=None, local=False, dyn=None):
        s = self.sim
        s0 = getattr(s, 'sea_level', None)
        saved_w = None
        if world is not None:
            saved_w = {k: s.world[k] for k in world}
            s.world.update(world)
        saved_f = (s.flex_dry, s.flex_wet)
        if local:
            s.flex_dry, s.flex_wet = LocalAiry(s.p['rho_mantle']), LocalAiry(s.p['rho_mantle']-s.p['rho_water'])
        saved_d = s.dyn
        if dyn is not None:
            s.dyn = dyn
        try:
            z = s.elevation(thick, sed)
        finally:
            s.dyn = saved_d
            s.flex_dry, s.flex_wet = saved_f
            if saved_w is not None:
                s.world.update(saved_w)
            if s0 is None:
                if hasattr(s, 'sea_level'):
                    del s.sea_level
            else:
                s.sea_level = s0
        return np.asarray(z, np.float64)

    def keep(self, name, z, thick=None, sed=None, cont=None):
        w = self.sim.world
        self.surfaces[name] = (z, np.asarray(w['thick'] if thick is None else thick, np.float64).copy(),
                               np.asarray(w['sed'] if sed is None else sed, np.float64).copy(),
                               np.asarray(w['cont'] if cont is None else cont, bool).copy())

    def _wrap(self):
        s = self
        sim = self.sim
        orig_world_from, orig_erode, orig_update = sim._world_from, sim._erode, sim._update_rasters

        def world_from(idx, claimed, record=True):
            out = orig_world_from(idx, claimed, record)
            if record and s.prev is not None:
                s.keep('composed_olddyn', s.zdiag(dyn=s.prev['dyn']))
                s.keep('composed', s.zdiag())
            return out

        def erode(th):
            w = sim.world
            s.keep('tectonics', s.zdiag(th, w['sed']), thick=th)
            e_cont, dsed = orig_erode(th)
            s.keep('eroded', s.zdiag(th+e_cont, w['sed']+dsed), thick=th+e_cont, sed=w['sed']+dsed)
            return e_cont, dsed

        def update(dthick, idx, dsed, drock):
            s._transfer_read(dthick, idx, dsed)
            out = orig_update(dthick, idx, dsed, drock)
            s.keep('transferred', s.zdiag())
            s._returned(idx)
            return out

        sim._world_from, sim._erode, sim._update_rasters = world_from, erode, update

    def unwrap(self):
        for k in ('_world_from', '_erode', '_update_rasters'):
            self.sim.__dict__.pop(k, None)

    def _transfer_read(self, dthick, idx, dsed):
        """Continental plate nodes no world node reads this step, and the increments they miss."""
        sim, g = self.sim, self.sim.g
        a = g.area
        new = sim.new_crust
        unread_a = cont_a = 0.
        miss_t, miss_s, app_t, app_s, wm, wa = [], [], [], [], [], []
        for k, pl in enumerate(sim.plates):
            ws = np.flatnonzero((sim.owner_slot == k) & ~new)
            js = idx[k, ws]
            read = np.bincount(js, minlength=g.size) > 0
            cn = pl.mask & pl.cont
            un = cn & ~read
            cont_a += a[cn].sum(); unread_a += a[un].sum()
            if un.any():
                wi = fields.nearest_nodes(g, g.xyz[un]@pl.Q.T)       # world position of plate node j: Q p_j
                ok = sim.owner_slot[wi] == k
                miss_t.append(dthick[wi][ok]); wm.append(a[un][ok])
                miss_s.append((dsed[wi] if dsed is not None else np.zeros(wi.size))[ok])
            rc = cn & read
            if rc.any():
                d_int = np.bincount(js, weights=dthick[ws], minlength=g.size)/np.maximum(np.bincount(js, minlength=g.size), 1)
                app_t.append(d_int[rc]); wa.append(a[rc])
                if dsed is not None:
                    d_s = np.bincount(js, weights=dsed[ws], minlength=g.size)/np.maximum(np.bincount(js, minlength=g.size), 1)
                    app_s.append(d_s[rc])
        rms = lambda v, w: float(np.sqrt(np.average(np.concatenate(v)**2, weights=np.concatenate(w)))) if v and np.concatenate(w).sum() > 0 else 0.
        self.transfer = {'unread_cont_share': float(unread_a/max(cont_a, 1.)),
                         'missed_rms_thick_m': rms(miss_t, wm), 'missed_rms_sed_m': rms(miss_s, wm),
                         'applied_rms_thick_m': rms(app_t, wa), 'applied_rms_sed_m': rms(app_s, wa) if app_s else 0.}

    def _returned(self, idx):
        sim, w = self.sim, self.sim.world
        th, sd, cn = w['thick'].copy(), w['sed'].copy(), w['cont'].copy()
        new = sim.new_crust
        live = np.zeros(sim.g.size, bool)
        for k, pl in enumerate(sim.plates):
            sel = np.flatnonzero((sim.owner_slot == k) & ~new)
            j = idx[k, sel]
            ok = pl.mask[j]
            sel, j = sel[ok], j[ok]
            th[sel] = pl.thick[j]; sd[sel] = pl.sed[j]; cn[sel] = pl.cont[j]; live[sel] = True
        z = self.zdiag(world={'thick': th, 'sed': sd, 'cont': cn})
        self.keep('returned', z, thick=th, sed=sd, cont=cn)
        self.transfer['returned_live_share'] = float(live.mean())

    def quiet_state(self):
        """state() without moving the trajectory: state() re-solves the sea, and the next solve's loads read it."""
        s = self.sim
        s0 = getattr(s, 'sea_level', None)
        st = s.state()
        if s0 is not None:
            s.sea_level = s0
        return st

    def marks(self, st):
        """Parcel marks at prev_end, before any selection on height."""
        sim, g, geo = self.sim, self.sim.g, self.geo
        w = sim.world
        T = w['thick']
        tcls = np.digitize(T, T_BINS[1:-1])
        de = geo.dist_to(~w['cont'], 700e3)
        ecls = np.digitize(de, E_BINS[1:-1])
        bt = st['boundary_type'] > 0
        self.bhist = (self.bhist+[bt])[-10:]
        near = geo.dist_to(np.any(self.bhist, axis=0), 151e3) <= 150e3
        return {'T': tcls, 'edge': ecls, 'boundary': near.astype(int)}

    def step(self):
        """One traced step; returns its record."""
        sim = self.sim
        st0 = self.prev_state
        self.surfaces = {}
        self.prev = {'dyn': sim.dyn.copy(), 'sea': float(getattr(sim, 'sea_level', 0.))}
        self.keep('prev_end', np.asarray(st0['elevation_m'], np.float64))
        mk = self.marks(st0)
        sim.step()
        st = sim.state()
        z_end = np.asarray(st['elevation_m'], np.float64)
        self.keep('end', z_end)
        self.keep('end_local', self.zdiag(local=True))
        self.keep('end_prev_sea', z_end+float(getattr(sim, 'sea_level', 0.))-self.prev['sea'])
        rec = {'t': float(sim.t), 'transfer': dict(self.transfer), 'stages': {}}
        prev_masks = None
        for name in STAGES+SIDE:
            if name not in self.surfaces:
                continue
            z, T, sd, cn = self.surfaces[name]
            r, m = self.geo.measures(z, T, sd, cn)
            base = prev_masks if name in STAGES else self._chain_masks.get(
                {'returned': 'transferred', 'end_local': 'end', 'end_prev_sea': 'end'}[name])
            if base is not None:
                a = sim.g.area
                for key in ('hole', 'island'):
                    made = m[key] & ~base[key]
                    r[f'{key}s_made'] = int(made.sum()); r[f'{key}s_gone'] = int((base[key] & ~m[key]).sum())
                    if key == 'hole':
                        r['holes_made_by_mark'] = {mn: np.bincount(mv[made], minlength=4 if mn != 'boundary' else 2).tolist()
                                                   for mn, mv in mk.items()}
                r['band_in'] = float(a[m['band'] & ~base['band']].sum()/self.geo.A)
                r['band_out'] = float(a[base['band'] & ~m['band']].sum()/self.geo.A)
                r['far_low_in'] = float(a[m['far_low'] & ~base['far_low']].sum()/self.geo.A)
                r['far_low_out'] = float(a[base['far_low'] & ~m['far_low']].sum()/self.geo.A)
            if name in STAGES:
                prev_masks = m
                self._chain_masks[name] = m
            rec['stages'][name] = r
        rec['mark_area_share'] = {mn: [float(sim.g.area[mv == c].sum()/self.geo.A) for c in range(4 if mn != 'boundary' else 2)]
                                  for mn, mv in mk.items()}
        self.prev_state = st
        return rec

    def run(self, steps, out, mode):
        self._chain_masks = {}
        self.prev_state = self.quiet_state()
        with open(out, 'w') as fh:
            for i in range(steps):
                t0 = time.perf_counter()
                self._chain_masks = {}
                rec = self.step()
                rec['mode'] = mode; rec['seconds'] = time.perf_counter()-t0
                fh.write(json.dumps(rec)+'\n'); fh.flush()
                e = rec['stages']['end']
                print(f"{mode} t={rec['t']:+.0f}: holes {e['holes']}, islands {e['islands']}, crossings "
                      f"{e['crossings_per_mkm2']:.1f}, lowflat4 {e['lowflat4_of_land']:.3f}, unread "
                      f"{rec['transfer']['unread_cont_share']:.3f} ({rec['seconds']:.0f} s)", flush=True)


def smooth_material(sim, km=SMOOTH_KM):
    """Conservative graph diffusion of every plate's continental crust and sediment over its continental nodes."""
    g = sim.g
    total = (km*1e3)**2/2
    steps = max(1, int(np.ceil(total*g.diffusion_bound()/.9)))
    for pl in sim.plates:
        dom = pl.mask & pl.cont
        if dom.sum() < 2:
            continue
        th = pl.thick.astype(np.float64); sd = pl.sed.astype(np.float64)
        for _ in range(steps):
            th = gl.diffuse(th, dom, total/steps, 1., g.nbr, g.dist, g.axis, g.area)
            sd = gl.diffuse(sd, dom, total/steps, 1., g.nbr, g.dist, g.axis, g.area)
        pl.thick[dom] = th[dom].astype(np.float32); pl.sed[dom] = np.maximum(sd[dom], 0.).astype(np.float32)
        pl.rock[dom] = np.minimum(pl.rock[dom], pl.sed[dom])
    sim.compose_from_plates()


def null_run(sim, steps, out, mode):
    """Rigid rotation and composition only: no dynamic change, no processes, no transfer."""
    geo = Geo(sim.g)
    w = sim.world
    dt = sim.p['dt_myr']
    st = sim.state()
    with open(out, 'w') as fh:
        z = np.asarray(st['elevation_m'], np.float64)
        prev = None
        for i in range(steps+1):
            if i:
                for pl in sim.plates:
                    pl.Q = fields.rotation(pl.axis, pl.speed*1e-2*dt*1e6/sim.g.radius)@pl.Q
                sim.t += dt
                idx, claimed = sim._claims()
                sim._world_from(idx, claimed, record=False)
                w = sim.world
                z = np.asarray(sim.elevation(), np.float64)
            r, m = geo.measures(z, w['thick'], w['sed'], w['cont'])
            if prev is not None:
                r['holes_made'] = int((m['hole'] & ~prev['hole']).sum()); r['holes_gone'] = int((prev['hole'] & ~m['hole']).sum())
            prev = m
            fh.write(json.dumps({'mode': mode, 't': float(sim.t), 'step': i, 'stages': {'end': r}})+'\n'); fh.flush()
            print(f"{mode} step {i}: holes {r['holes']}, islands {r['islands']}, crossings {r['crossings_per_mkm2']:.1f}, "
                  f"rough_T {r['rough_T']:.0f}, rough_sed {r['rough_sed']:.0f}", flush=True)


def warm(out, seed, stop, n, years, params):
    out.mkdir(parents=True, exist_ok=True)
    g = CubeSphere(n)
    sim = PlateSim(g, seed=seed, t0=-years, params=params)
    sim.state()
    t0 = time.perf_counter()
    while sim.t < stop-1e-6:
        sim.step()
        sim.state()
        if round(sim.t) % 10 == 0:
            print(f'warm seed {seed} t={sim.t:+.0f} ({time.perf_counter()-t0:.0f} s)', flush=True)
    pickle_save(sim, out/'warm.pkl')
    (out/'warm.json').write_text(json.dumps({'seed': seed, 'n': n, 'years': years, 'stop': stop, 'params': params,
                                             't': sim.t, 'seconds': time.perf_counter()-t0}, indent=1))


def summary(out):
    """Per mode: stage means over the traced steps (the reading table for journal 142's rule)."""
    keys = ('land', 'holes', 'islands', 'crossings_per_mkm2', 'dry4_of_land', 'lowflat4_of_land', 'low200_of_far300',
            'far300_of_land', 'band_of_globe', 'rough_z', 'rough_T', 'rough_sed', 'holes_made', 'holes_gone',
            'islands_made', 'islands_gone', 'band_in', 'band_out', 'far_low_in', 'far_low_out')
    lines = []
    for f in sorted(Path(out).glob('trace-*.jsonl')):
        recs = [json.loads(l) for l in f.read_text().splitlines() if l.strip()]
        if not recs:
            continue
        mode = recs[0]['mode']
        lines.append(f'\n## {mode} ({len(recs)} steps, t {recs[0]["t"]:+.0f} ... {recs[-1]["t"]:+.0f})')
        names = [n for n in STAGES+SIDE if n in recs[-1]['stages']]
        lines.append('stage | ' + ' | '.join(keys))
        for n in names:
            rs = [r['stages'][n] for r in recs if n in r['stages']]
            vals = []
            for k in keys:
                v = [r[k] for r in rs if k in r]
                vals.append(f'{np.mean(v):.4g}' if v else '-')
            lines.append(f'{n} | ' + ' | '.join(vals))
        if 'transfer' in recs[0] and recs[0]['transfer']:
            tk = recs[0]['transfer'].keys()
            lines.append('transfer: ' + ', '.join(f"{k} {np.mean([r['transfer'][k] for r in recs if k in r['transfer']]):.4g}" for k in tk))
        made = {}
        for r in recs:
            for n, s in r['stages'].items():
                for mn, v in s.get('holes_made_by_mark', {}).items():
                    made.setdefault((n, mn), np.zeros(len(v)))
                    made[(n, mn)] += v
        if made:
            lines.append('holes made by mark (sum over steps): stage | T <30/30-38/38-45/>=45 km | edge <100/100-300/300-600/>=600 km | boundary no/yes')
            for n in names:
                if (n, 'T') in made:
                    lines.append(f"{n} | {made[(n, 'T')].astype(int).tolist()} | {made[(n, 'edge')].astype(int).tolist()} | "
                                 f"{made[(n, 'boundary')].astype(int).tolist()}")
            sh = recs[0].get('mark_area_share')
            if sh:
                lines.append(f"mark area shares (globe, first step): {json.dumps(sh)}")
    text = '\n'.join(lines)
    (Path(out)/'summary.txt').write_text(text)
    print(text)


def check_identity(out, n, steps, seed):
    """A traced run equals an untraced one, field by field (threads 1)."""
    params = {'threads': 1}
    g = CubeSphere(n)
    a = PlateSim(g, seed=seed, t0=-20., params=params)
    b = PlateSim(g, seed=seed, t0=-20., params=params)
    a.state(); b.state()
    for _ in range(3):                        # a few steps first (the first steps carry no world_from record)
        a.step(); a.state(); b.step(); b.state()
    tr = Tracer(b)
    tr.prev_state = tr.quiet_state()
    tr._chain_masks = {}
    for i in range(steps):
        a.step(); sa = a.state()
        tr._chain_masks = {}
        tr.step(); sb = tr.prev_state
        bad = [k for k in sa if not np.array_equal(np.asarray(sa[k]), np.asarray(sb[k]), equal_nan=True)]
        print(f'step {i}: {"IDENTICAL" if not bad else "DIFFER " + ",".join(bad)}', flush=True)
        if bad:
            return False
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=('warm', 'trace', 'summary', 'check-identity'))
    ap.add_argument('out')
    ap.add_argument('--seed', type=int, default=1028)
    ap.add_argument('--stop', type=float, default=-20.)
    ap.add_argument('--years', type=float, default=250.)
    ap.add_argument('--n', type=int, default=217)
    ap.add_argument('--mode', default='default', choices=('default', 'smoothed', 'null', 'null_raw'))
    ap.add_argument('--steps', type=int, default=20)
    ap.add_argument('--param', action='append', default=[], help='key=json value (e.g. threads=1)')
    a = ap.parse_args(argv)
    out = Path(a.out)
    params = {'threads': 1}
    for kv in a.param:
        k, v = kv.split('=', 1)
        params[k] = json.loads(v)
    if a.cmd == 'warm':
        warm(out, a.seed, a.stop, a.n, a.years, params)
    elif a.cmd == 'trace':
        sim = pickle_load(out/'warm.pkl')
        f = out/f'trace-{a.mode}.jsonl'
        if a.mode in ('smoothed', 'null'):
            smooth_material(sim)
        if a.mode in ('null', 'null_raw'):
            null_run(sim, a.steps, f, a.mode)
        else:
            Tracer(sim).run(a.steps, f, a.mode)
    elif a.cmd == 'summary':
        summary(out)
    else:
        sys.exit(0 if check_identity(out, a.n, a.steps, a.seed) else 1)


if __name__ == '__main__':
    main()
