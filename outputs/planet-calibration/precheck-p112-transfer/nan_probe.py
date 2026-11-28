"""Pass 112: where the repaired C goes non-finite (Transfer Sweep, seed 1004, between -140 and -130 Myr), calling state() after every step as simulate_plates does.
Replays the seed with C, checks every `_tent_return` call's input and output and the plate fields after each step;
at the first non-finite value prints the call and saves its inputs to nan_case.npz. Run from the repo root:
  .venv\\Scripts\\python.exe -B outputs\\planet-calibration\\precheck-p112-transfer\\nan_probe.py --seed 1004
"""
import argparse
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=1004)
    ap.add_argument('--steps', type=int, default=130)
    a = ap.parse_args()
    g = CubeSphere(217)
    sim = PlateSim(g, seed=a.seed, t0=-250., params={'remap_rule': 'overlap', 'rift_area': 'plates', 'threads': 2})   # as the sweep (NUMBA_NUM_THREADS=2 too)
    arch = Path('outputs/planet-calibration/sweep-20261127-173125/remap_rule-overlap_rift_area-plates')/f'seed{a.seed}'/'snapshots'
    sim.state()
    orig = sim._tent_return

    def checked(k, pl, wn, own, d, name):
        cur = sim._return_cur(pl, name)
        J, out = orig(k, pl, wn, own, d, name)
        ws, cs, f = sim._tent[k]
        bad_in = ~np.isfinite(d[ws]).all() or ~np.isfinite(cur[pl.mask]).all()
        if bad_in or not np.isfinite(out).all():
            print(f't {sim.t} plate slot {k} field {name}: input finite {not bad_in}, output non-finite '
                  f'{int((~np.isfinite(out)).sum())} of {out.size}', flush=True)
            np.savez(HERE/'nan_case.npz', d=d, ws=ws, cs=cs, f=f, cur=cur.astype(np.float64), mask=pl.mask,
                     cont=pl.cont, wn=wn, own=own, new_crust=sim.new_crust, wcont=sim.world['cont'],
                     wtype=sim._tent_type[k], J=J, out=out, area=g.area)
            raise SystemExit(1)
        return J, out
    sim._tent_return = checked
    for s in range(a.steps):
        sim.step()
        for pl in sim.plates:
            for nm in ('thick', 'sed', 'rock', 'therm'):
                v = getattr(pl, nm)[pl.mask]
                if nm != 'therm' and not np.isfinite(v).all():
                    print(f'after step t {sim.t}: plate {pl.id} field {nm} non-finite', flush=True)
                    raise SystemExit(1)
        st = sim.state()                             # as simulate_plates' record (the sea level feeds the next step)
        if s % 10 == 9:
            f = arch/f't{sim.t:+07.1f}.npz'
            same = None
            if f.exists():
                with np.load(f) as z:
                    same = all(np.array_equal(z[k], v, equal_nan=True) for k, v in st.items() if k in z.files)
            print('t', sim.t, 'ok; matches the sweep snapshot:', same, flush=True)
    print('no non-finite value')


if __name__ == '__main__':
    main()
