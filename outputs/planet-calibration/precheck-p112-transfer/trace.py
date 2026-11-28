"""Pass 112 (journal 153), rule R2 / R3: C (overlap + plates) with the material-aware return, traced per step.

Adapted from the review's remap_trace.py (outputs/review-astra/2026-11-27/): the observer wraps `_tent_return` and
returns its arrays unchanged. Per step it writes the oceanic nodes that cross below 6,999 m in the thickness return,
the thinnest oceanic node of the state, the oceanic deficit below 7 km (m^3), the requested thickness loss, and the
booked shortfall / unplaced / clamped volumes. Run from the repo root:

  .venv\\Scripts\\python.exe -B outputs\\planet-calibration\\precheck-p112-transfer\\trace.py --seed 1001 --steps 10
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim

HERE = Path(__file__).resolve().parent
KEYS = ('return_shortfall_thick_m3', 'return_unplaced_thick_m3', 'return_thick_clamped_m3', 'return_rescaled_thick_m3',
        'return_shortfall_sed_m3', 'return_sed_clamped_m3')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=10)
    ap.add_argument('--seed', type=int, default=1001)
    ap.add_argument('--out', type=Path, default=None)
    args = ap.parse_args()
    out_path = args.out or HERE/f'trace_{args.seed}_{args.steps}.csv'
    g = CubeSphere(217)
    sim = PlateSim(g, seed=args.seed, params={'remap_rule': 'overlap', 'rift_area': 'plates', 'threads': 2})
    floor = sim.p['ocean_thick_m']
    original = sim._tent_return
    step = {}

    def observed(k, pl, wn, own, d, name):
        J, out = original(k, pl, wn, own, d, name)
        if name == 'thick':
            ws = sim._tent[k][0]
            step['req_loss'] += float((g.area[ws]*np.minimum(d[ws], 0.)).sum())
            oc = ~pl.cont[J]
            before = pl.thick[J].astype(np.float64)
            new = oc & (before >= floor-1) & (before+out < floor-1)
            step['cross'] += int(new.sum())
            step['cross_km3'] += float((np.maximum(floor-1-before[new]-out[new], 0)*g.area[J[new]]).sum())/1e9
        return J, out
    sim._tent_return = observed
    sim.state()
    rows = []
    for s in range(args.steps):
        step.update(req_loss=0., cross=0, cross_km3=0.)
        ET0 = {k: sim.erosion_totals.get(k, 0.) for k in KEYS}
        sim.step()
        st = sim.state()
        oc = ~st['continental']
        th = st['crust_thickness_m'][oc]
        row = dict(time=sim.t, cross=step['cross'], cross_km3=step['cross_km3'],
                   oc_min_m=float(th.min()),
                   oc_deficit_km3=float((np.clip(floor-th, 0, None)*g.area[oc]).sum())/1e9,
                   req_loss_km3=step['req_loss']/1e9)
        for k in KEYS:
            row[k.replace('_m3', '_km3')] = (sim.erosion_totals.get(k, 0.)-ET0[k])/1e9
        rows.append(row)
        if s % 10 == 9 or s == args.steps-1:
            print(s+1, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    with out_path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)
    tot = {k: sum(r[k] for r in rows) for k in rows[0] if k not in ('time', 'oc_min_m', 'oc_deficit_km3')}
    print('totals', tot)
    print('end: oceanic min', rows[-1]['oc_min_m'], 'deficit km3', rows[-1]['oc_deficit_km3'])
    print('shortfall / requested loss (thick):', tot['return_shortfall_thick_km3']/max(-tot['req_loss_km3'], 1e-30))


if __name__ == '__main__':
    main()
