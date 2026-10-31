"""Paired probe for force_law 'torque' (forty-third pass): the same plate states continued under both force laws.

A control run (PARAMS, threads 1) is stepped to each start time; there the simulator is copied twice: A goes on
unchanged, B switches to force_law 'torque' (its plates start from A's motions and relax toward the torque balance
over tau_speed_myr). Both run --span Myr; every --every Myr the contour lengths and rates
(planet_reviews.trench_origin.snapshot_metrics) are printed: intra-oceanic and continental-side trench, ridge
production, the closing p50 of intra-oceanic trench and the opening p50 of ridges, and the NNR speeds of plates over
10 Mkm^2 with little slab (slab fraction of the perimeter < 0.05) and with much (> 0.25); then the paired differences
B - A at the end and as span means.

  python tools/torque_paired_probe.py --seed 1001 --starts=-200,-150 --span 30 --every 5 [--gate] --out <folder>

--gate adds arm C: torque + subduction_initiation 'convergence' (pass 42's slab gate), every trench of the start given
its slab, so only convergence that begins after the start has to build one; printed with C - B.
"""
import argparse
import copy
import copyreg
import importlib
import json
import types
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_reviews import trench_origin as to
from planet_tectonics.plates import PlateSim, SUBDUCTION
from planet_tectonics.motion import SLAB_BAND_RINGS
from planet_grid import fields

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))   # the flexure keeps a module


def seed_slabs(sim):
    """Every overriding subduction node of this moment gets its slab (as tools/slab_paired_probe.py)."""
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


def metrics(sim, tmp, t0=None):
    save(sim, tmp)
    r = to.snapshot_metrics(sim.g, tmp)
    big = [q for q in r['plates_rows'] if q['area_mkm2'] > 10]
    def wmean(qs):
        return float(np.average([q['nnr_cm_yr'] for q in qs], weights=[q['area_mkm2'] for q in qs])) if qs else float('nan')
    return {k: r[k] for k in ('io_km', 'cont_sub_km', 'production_km2_yr', 'ridge_km', 'plates', 'ocb_km')} | {
        'ocb_p50': r['ocb_rate_p10_p50_p90'][1],
        'io_p50': r['io_rate_p10_p50_p90'][1], 'ridge_p50': r['ridge_rate_p50'],
        'v_slabless': wmean([q for q in big if q['slab_fraction'] < .05]),
        'v_slab': wmean([q for q in big if q['slab_fraction'] > .25]),
        'v_all': wmean(big)} | {f'ev_{k}': sum(1 for e in sim.events if e['t'] > t0 and e['kind'] == k)
                                for k in ('rift', 'subduction initiation', 'detach', 'merge')}


KEYS = (('io_km', 1e3, 'k km'), ('cont_sub_km', 1e3, 'k km'), ('production_km2_yr', 1., 'km^2/yr'), ('ridge_km', 1e3, 'k km'),
        ('io_p50', 1., 'mm/yr'), ('ridge_p50', 1., 'mm/yr'), ('v_slabless', 1., 'cm/yr'), ('v_slab', 1., 'cm/yr'),
        ('v_all', 1., 'cm/yr'), ('plates', 1., ''), ('ev_detach', 1., 'events'), ('ev_merge', 1., 'events'))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--seed', type=int, default=1001)
    ap.add_argument('--n', type=int, default=217)
    ap.add_argument('--starts', default='-200,-150,-100')
    ap.add_argument('--span', type=int, default=30)
    ap.add_argument('--every', type=int, default=5)
    ap.add_argument('--out', required=True)
    ap.add_argument('--gate', action='store_true', help='add arm C: torque + subduction_initiation convergence, every '
                    'trench of the start given its slab (as tools/slab_paired_probe.py)')
    a = ap.parse_args(argv)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tmp = out/'tmp'
    g = CubeSphere(a.n)
    base = PlateSim(g, seed=a.seed, params={'threads': 1})
    res = []
    for t0 in sorted(float(s) for s in a.starts.split(',')):
        while base.t < t0-1e-9:
            base.step()
        A, B = copy.deepcopy(base), copy.deepcopy(base)
        B.p['force_law'] = 'torque'
        arms = {'A': A, 'B': B}
        if a.gate:
            C = copy.deepcopy(base)
            C.p['force_law'] = 'torque'; C.p['subduction_initiation'] = 'convergence'
            seed_slabs(C)
            arms['C'] = C
        rows = []
        for k in range(1, a.span+1):
            for sim in arms.values():
                sim.step()
            if k % a.every == 0:
                M = {x: metrics(sim, tmp, t0) for x, sim in arms.items()}
                rows.append({'t': A.t, **M})
                f2 = lambda key, sc=1e3, fmt='{:5.1f}': ' / '.join(fmt.format(M[x][key]/sc) for x in M)
                print(f"start {t0:.0f}  t {A.t:6.0f}  io {f2('io_km')}k ({f2('io_p50', 1., '{:.0f}')} mm/yr)  cs {f2('cont_sub_km')}k"
                      f"  prod {f2('production_km2_yr', 1., '{:.2f}')}  ridge {f2('ridge_km')}k ({f2('ridge_p50', 1., '{:.0f}')})"
                      f"  v slab-less {f2('v_slabless', 1., '{:.1f}')} slab {f2('v_slab', 1., '{:.1f}')} all {f2('v_all', 1., '{:.1f}')}"
                      f"  plates {f2('plates', 1., '{:.0f}')}  events rift/init/detach/merge "
                      f"{' | '.join('/'.join(str(M[x][f'ev_{k}']) for k in ('rift', 'subduction initiation', 'detach', 'merge')) for x in M)}"
                      f"  OCB {f2('ocb_km')}k (p50 {f2('ocb_p50', 1., '{:.0f}')})", flush=True)
        res.append({'start': t0, 'rows': rows})
        (out/'paired.json').write_text(json.dumps(res, indent=1))
    pairs = (('B', 'A'),)+((('C', 'B'),) if a.gate else ())
    print('paired differences B - A (torque - lengths)' + (', C - B (the slab gate under torque)' if a.gate else '') +
          ': end of each span; span mean')
    for key, sc, u in KEYS+(('ocb_km', 1e3, 'k km'),):
        for x, y in pairs:
            end = [(r['rows'][-1][x][key]-r['rows'][-1][y][key])/sc for r in res]
            mean = [np.nanmean([q[x][key]-q[y][key] for q in r['rows']])/sc for r in res]
            print(f"  {key} {x}-{y}: end {' / '.join(f'{v:+.2f}' for v in end)}; span mean {' / '.join(f'{v:+.2f}' for v in mean)}"
                  f" (all: {np.nanmean(mean):+.2f}) {u}")


if __name__ == '__main__':
    main()
