"""Paired probe for subduction_initiation (forty-second pass): the same plate states continued under both rules.

A control run (PARAMS, threads 1) is stepped to each start time; there the simulator is copied three times: A goes on
unchanged, B switches to trench_polarity pair, C to pair + subduction_initiation convergence with every trench of
that moment given its slab (the band behind each overriding subduction node set to slab_init_km), so only convergence
that begins after the start time has to build one. All run --span Myr; every --every Myr the contour lengths
(planet_reviews.trench_origin.snapshot_metrics: intra-oceanic and continental-side trench with a slab, convergence
without a slab, ridge production) are printed, and the paired differences B - A (the polarity fix) and C - B (the slab
rule) at the end.

  python tools/slab_paired_probe.py --seed 1001 --starts -150,-100,-50 --span 30 --every 5 --out <folder>
"""
import argparse
import copy
import copyreg
import importlib
import json
import types
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import trench_origin as to
from planet_tectonics.motion import SLAB_BAND_RINGS
from planet_tectonics.plates import PlateSim, SUBDUCTION

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))   # the flexure keeps a module


def seed_slabs(sim):
    """Every overriding subduction node of this moment gets its slab: its band on the plate raster set to slab_init_km."""
    g, D = sim.g, sim.p['slab_init_km']
    over = (sim.btype == SUBDUCTION) & ~sim.trench
    band = sim._spread_same_plate(over.astype(float), SLAB_BAND_RINGS) > 0
    byid = {pl.id: pl for pl in sim.plates}
    for k in np.unique(sim.owner[band]):
        pl = byid.get(int(k))
        if pl is None:
            continue
        I = np.flatnonzero(band & (sim.owner == k))
        r = fields.nearest_nodes(g, g.xyz[I]@pl.Q)
        r = r[pl.mask[r]]
        pl.slab[r] = np.maximum(pl.slab[r], D)


def save(sim, f):
    st = sim.state()
    np.savez(str(f)+'.npz', **{k: np.asarray(v) for k, v in st.items()})
    Path(str(f)+'-plates.json').write_text(json.dumps(sim.plate_table()))


def metrics(sim, tmp):
    save(sim, tmp)
    r = to.snapshot_metrics(sim.g, tmp)
    return {k: r[k] for k in ('io_km', 'cont_sub_km', 'ocb_km', 'noslab_cont_km', 'production_km2_yr', 'ridge_km')} | {
        'plates': r['plates'], 'ocb_p50': r['ocb_rate_p10_p50_p90'][1]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--seed', type=int, default=1001)
    ap.add_argument('--n', type=int, default=217)
    ap.add_argument('--starts', default='-150,-100,-50')
    ap.add_argument('--span', type=int, default=30)
    ap.add_argument('--every', type=int, default=5)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tmp = out/'tmp'
    g = CubeSphere(a.n)
    base = PlateSim(g, seed=a.seed, params={'threads': 1})
    res = []
    for t0 in sorted(float(s) for s in a.starts.split(',')):
        while base.t < t0-1e-9:
            base.step()
        A, B, C = copy.deepcopy(base), copy.deepcopy(base), copy.deepcopy(base)
        B.p['trench_polarity'] = 'pair'
        C.p['trench_polarity'] = 'pair'
        C.p['subduction_initiation'] = 'convergence'
        seed_slabs(C)
        rows = []
        for k in range(1, a.span+1):
            A.step(); B.step(); C.step()
            if k % a.every == 0:
                M = {'A': metrics(A, tmp), 'B': metrics(B, tmp), 'C': metrics(C, tmp)}
                rows.append({'t': A.t, **M})
                f3 = lambda key, sc=1e3, fmt='{:5.1f}': ' / '.join(fmt.format(M[x][key]/sc) for x in 'ABC')
                print(f"start {t0:.0f}  t {A.t:6.0f}  io {f3('io_km')}k  cs {f3('cont_sub_km')}k  OCB {M['C']['ocb_km']/1e3:4.1f}k "
                      f"(p50 {M['C']['ocb_p50']:.0f})  prod {f3('production_km2_yr', 1., '{:.2f}')}  ridge {f3('ridge_km')}k  plates "
                      f"{f3('plates', 1., '{:.0f}')}", flush=True)
        res.append({'start': t0, 'rows': rows})
        (out/'paired.json').write_text(json.dumps(res, indent=1))
    print('paired differences (end of each span; mean over the span): B - A (pair polarity), C - B (slab rule)')
    for key, sc, u in (('io_km', 1e3, 'k km'), ('cont_sub_km', 1e3, 'k km'), ('production_km2_yr', 1., 'km^2/yr'),
                       ('ridge_km', 1e3, 'k km'), ('ocb_km', 1e3, 'k km')):
        for x, y in (('B', 'A'), ('C', 'B')):
            end = [(r['rows'][-1][x][key]-r['rows'][-1][y][key])/sc for r in res]
            mean = [np.mean([q[x][key]-q[y][key] for q in r['rows']])/sc for r in res]
            print(f"  {key} {x}-{y}: end {' / '.join(f'{v:+.2f}' for v in end)}; span mean {' / '.join(f'{v:+.2f}' for v in mean)}"
                  f" (all: {np.mean(mean):+.2f}) {u}")


if __name__ == '__main__':
    main()
