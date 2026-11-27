"""Pass 65: where node-scale crust spikes (the towers' elastic part) come from. Today's PARAMS (plus --param
overrides), one seed, every step: the highest node, its elevation split into Airy + dyn + plate support, its crust
against its neighbours' mean, and the thickening it got this step from each process (flow = the diffusion,
rift, subduction band, collision volume, erosion) plus the jump when the world re-reads the plate rasters
(_world_from: thickness after it minus the world's thickness at the end of the previous step).
Also per step: the largest spike anywhere (crust minus neighbours' mean) and what fed it."""
import sys, json, time, argparse
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim

ap = argparse.ArgumentParser()
ap.add_argument('--seed', type=int, default=1007)
ap.add_argument('--years', type=float, default=250)
ap.add_argument('--param', action='append', default=[])
ap.add_argument('--out', default='probe.jsonl')
a = ap.parse_args()
params = {'threads': 1}
for kv in a.param:
    k, v = kv.split('=')
    try:
        v = json.loads(v)
    except Exception:
        pass
    params[k] = v
g = CubeSphere(217)
sim = PlateSim(g, seed=a.seed, t0=-a.years, params=params)
rc, rm, rs = sim.p['rho_crust'], sim.p['rho_mantle'], sim.p['rho_sediment']
nbr = g.nbr
valid = nbr >= 0
rec = {}

def wrap(name, key, pick=lambda r: r):
    f = getattr(sim, name)
    def w(*args, **kw):
        r = f(*args, **kw)
        rec[key] = np.array(pick(r), dtype=float).copy()
        return r
    setattr(sim, name, w)

wrap('_rift', 'rift')
wrap('_subduction_thickening', 'sub')
wrap('_erode', 'erode', lambda r: r[0])
f_world = sim._world_from
def world_from(*args, **kw):
    prev = sim.world['thick'].copy() if hasattr(sim, 'world') and 'thick' in sim.world else None
    pc = sim.world['cont'].copy() if prev is not None else None
    r = f_world(*args, **kw)
    rec['map'] = np.where(pc & sim.world['cont'], sim.world['thick']-prev, 0.) if prev is not None else np.zeros(g.size)
    return r
sim._world_from = world_from
f_upd = sim._update_rasters
def upd(*args, **kw):
    pre = sim.tectonic.copy()
    r = f_upd(*args, **kw)
    rec['coll'] = sim.tectonic-pre
    rec['tect_pre'] = pre
    return r
sim._update_rasters = upd

def nbr_mean(x, cont=None):
    m = valid if cont is None else valid & cont[np.clip(nbr, 0, None)]
    xs = np.where(m, x[np.clip(nbr, 0, None)], 0.)
    return xs.sum(1)/np.maximum(m.sum(1), 1), m.sum(1)

out = open(a.out, 'w')
t0 = time.time()
steps = int(round(a.years))
for s in range(steps):
    rec.clear()
    sim.step()
    w = sim.world
    z = sim.elevation()
    T = w['thick'].astype(float)
    from planet_tectonics.surface import sediment_column
    sh = sediment_column(w['sed'], w['cont'], sim.p)[0]
    S = np.where(w['cont'], sh, 0.)
    airy = sim.p['cont_ref_elev_m']+(T-sim.p['cont_ref_thick_m'])*(rm-rc)/rm+S*(1-rs/rm)-sim.sea_level
    dyn = np.asarray(sim.dyn, float) if np.ndim(sim.dyn) else np.zeros(g.size)
    sup = z-airy-dyn
    cont = w['cont']
    Tn, nc = nbr_mean(T, cont)
    spike = np.where(cont & (nc >= 4), T-Tn, -1e9)
    flow = rec['tect_pre']-rec['rift']-rec['sub']
    comp = {'flow': flow, 'rift': rec['rift'], 'sub': rec['sub'], 'coll': rec['coll'], 'erode': rec['erode'], 'map': rec['map']}
    def at(i):
        return {k: round(float(v[i])) for k, v in comp.items()}
    i = int(np.argmax(np.where(cont, z, -1e9)))
    j = int(np.argmax(spike))
    row = {'t': round(sim.t, 1), 'top': {'i': i, 'z': round(float(z[i])), 'airy': round(float(airy[i])), 'dyn': round(float(dyn[i])),
           'plate': round(float(sup[i])), 'T': round(float(T[i])), 'Tn': round(float(Tn[i])), **at(i)},
           'spike': {'j': j, 'dT': round(float(spike[j])), 'T': round(float(T[j])), 'z': round(float(z[j])), **at(j)},
           'n_plate2km': int(((sup > 2000) & cont).sum()), 'n_spike10km': int((spike > 10000).sum()),
           'max_flow': round(float(flow.max())), 'max_coll': round(float(rec['coll'].max())), 'max_sub': round(float(rec['sub'].max())),
           'max_map': round(float(rec['map'].max())), 'max_rift': round(float(rec['rift'].max())),
           'sec': round(time.time()-t0)}
    out.write(json.dumps(row)+'\n'); out.flush()
    if s % 10 == 0:
        print(row['t'], row['top']['z'], row['top']['plate'], row['spike']['dT'], row['sec'], flush=True)
