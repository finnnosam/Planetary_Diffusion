"""P-h (journal 66, section 5): paired probe for basal_traction. The same plate states continued with and without
traction toward the slab flow.

A base run (dyn_rule 'slabs', force_law 'torque', threads 1) is stepped to each start time (late enough for the slab
record to have filled: the first ~60 Myr are its transient); there the simulator is copied twice: A goes on unchanged,
B switches on basal_traction. Both run --span Myr; every --every Myr it prints C(d) of different plates (NNR, the
0-500, 500-1000 and 1000-2000 km bins and their mean), the intra-oceanic trench (T2's contour length) and the plate
count; then the paired differences B - A as span means.

Rule (fixed 31 Oct): sweep traction only if the span-mean C(d) < 2,000 km rises by >= +0.06 (half the Earth slab
proxy's +0.12) without T2 rising by > 2k km or the plate count jumping.

  python tools/traction_paired_probe.py --seed 1001 --starts=-120,-60 --span 30 --every 5 --param slab_alpha=depth --out <folder>
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
from planet_reviews import force_balance_earth as fbe
from planet_reviews import trench_origin as to
from planet_reviews import trench_pairs as tp
from planet_tectonics.plates import PlateSim

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))   # the flexure keeps a module


def save(sim, f):
    st = sim.state()
    np.savez(str(f)+'.npz', **{k: np.asarray(v) for k, v in st.items()})
    Path(str(f)+'-plates.json').write_text(json.dumps(sim.plate_table()))


def metrics(sim, tmp):
    save(sim, tmp)
    r = to.snapshot_metrics(sim.g, tmp)
    g = sim.g
    om = {pl.id: pl.axis*pl.speed*1e-2/g.radius for pl in sim.plates}
    w = fbe.remove_net_rotation(g, sim.owner, om)
    wn = np.stack([w.get(int(k), np.zeros(3)) for k in sim.owner])
    coh = tp.coherence(g, sim.owner, np.cross(wn, g.xyz), np.random.default_rng(0))
    v = sim.mantle.flow.v_xyz
    return {'coh': coh, 'coh2000': float(np.mean(coh[:3])), 'io_km': r['io_km'], 'cont_sub_km': r['cont_sub_km'],
            'plates': r['plates'], 'flow_rms_cm_yr': float(np.sqrt(np.average((v**2).sum(1), weights=g.area))*100)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--seed', type=int, default=1001)
    ap.add_argument('--n', type=int, default=217)
    ap.add_argument('--starts', default='-120,-60')
    ap.add_argument('--span', type=int, default=30)
    ap.add_argument('--every', type=int, default=5)
    ap.add_argument('--out', required=True)
    ap.add_argument('--param', action='append', default=[], metavar='KEY=VALUE', help='more base parameters (strings), '
                    'e.g. --param slab_alpha=depth')
    a = ap.parse_args(argv)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tmp = out/'tmp'
    g = CubeSphere(a.n)
    extra = dict(kv.split('=', 1) for kv in a.param)
    base = PlateSim(g, seed=a.seed, params={'threads': 1, 'dyn_rule': 'slabs', 'force_law': 'torque', **extra})
    res = []
    for t0 in sorted(float(s) for s in a.starts.split(',')):
        while base.t < t0-1e-9:
            base.step()
        A, B = copy.deepcopy(base), copy.deepcopy(base)
        B.p['basal_traction'] = True
        M0 = {x: metrics(s, tmp) for x, s in (('A', A), ('B', B))}
        print(f"start {t0:.0f}: C(d) <2,000 km {M0['A']['coh2000']:+.2f}, io {M0['A']['io_km']/1e3:.1f}k km, plates "
              f"{M0['A']['plates']}, flow rms {M0['A']['flow_rms_cm_yr']:.2f} cm/yr", flush=True)
        rows = []
        for k in range(1, a.span+1):
            A.step(); B.step()
            if k % a.every == 0:
                M = {'A': metrics(A, tmp), 'B': metrics(B, tmp)}
                rows.append({'t': A.t, **M})
                print(f"start {t0:.0f}  t {A.t:6.0f}  C(d) 0-500/500-1000/1000-2000 A "
                      + '/'.join(f'{c:+.2f}' for c in M['A']['coh'][:3])+'  B '+'/'.join(f'{c:+.2f}' for c in M['B']['coh'][:3])
                      + f"  mean {M['A']['coh2000']:+.2f} / {M['B']['coh2000']:+.2f}  io {M['A']['io_km']/1e3:.1f} / "
                      f"{M['B']['io_km']/1e3:.1f}k km  cs {M['A']['cont_sub_km']/1e3:.1f} / {M['B']['cont_sub_km']/1e3:.1f}k  "
                      f"plates {M['A']['plates']} / {M['B']['plates']}  flow {M['A']['flow_rms_cm_yr']:.2f} / "
                      f"{M['B']['flow_rms_cm_yr']:.2f} cm/yr", flush=True)
        res.append({'start': t0, 'start_metrics': M0['A'], 'rows': rows})
        (out/'paired.json').write_text(json.dumps(res, indent=1))
    print('\npaired B - A (traction - no traction), span means:')
    allc, alli, allp = [], [], []
    for r in res:
        dc = np.mean([q['B']['coh2000']-q['A']['coh2000'] for q in r['rows']])
        di = np.mean([q['B']['io_km']-q['A']['io_km'] for q in r['rows']])/1e3
        dp = np.mean([q['B']['plates']-q['A']['plates'] for q in r['rows']])
        allc.append(dc); alli.append(di); allp.append(dp)
        print(f"  start {r['start']:.0f}: C(d) <2,000 km {dc:+.3f}; io {di:+.1f}k km; plates {dp:+.1f}")
    ok = np.mean(allc) >= .06 and np.mean(alli) <= 2.
    print(f"  all starts: C(d) {np.mean(allc):+.3f} (rule >= +0.06), io {np.mean(alli):+.1f}k km (rule <= +2k), plates "
          f"{np.mean(allp):+.1f}  ->  P-h {'passes: sweep a traction arm' if ok else 'fails: traction is not swept'}")


if __name__ == '__main__':
    main()
