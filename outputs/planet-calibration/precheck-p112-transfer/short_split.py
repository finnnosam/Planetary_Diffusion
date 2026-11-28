"""Pass 112: where the booked thickness shortfall of the live-renormalised return comes from (R3 failed: 0.15% of the
returned loss on seed 1001 over 250 Myr). Replays K steps of C (as the sweep: threads 2, state() every step) and splits
each thickness return's requested loss into: cells with no live contributor (every contributor hidden or over new
crust), cells whose live contributors lack the material, and the cap after the rescaling; also the gain on cells with
no live contributor, and the same cells' contributors split into hidden (not owned) and over new crust.
  .venv\\Scripts\\python.exe -B outputs\\planet-calibration\\precheck-p112-transfer\\short_split.py --seed 1001 --steps 60
"""
import argparse

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=1001)
    ap.add_argument('--steps', type=int, default=60)
    a = ap.parse_args()
    g = CubeSphere(217)
    sim = PlateSim(g, seed=a.seed, t0=-250., params={'remap_rule': 'overlap', 'rift_area': 'plates', 'threads': 2})
    sim.state()
    orig = sim._tent_return
    acc = {}

    def split(k, pl, wn, own, d, name):
        if name == 'thick':
            ws, cs, f = sim._tent[k]
            live = pl.mask & own & ~sim.new_crust[wn]
            fl = f*live[cs]
            sf = fl.sum(1)
            has = sf > 0
            fl = fl/np.where(has, sf, 1.)[:, None]
            avail = np.maximum(pl.thick.astype(np.float64)-sim._return_floor(pl, 'thick'), 0.)
            ebar = (fl*avail[cs]).sum(1)
            dw = d[ws]
            loss, gain = np.minimum(dw, 0.), np.maximum(dw, 0.)
            aw = g.area[ws]
            hidden = (pl.mask & ~own)[cs].all(1)
            acc['loss'] = acc.get('loss', 0.)-float((aw*loss).sum())
            acc['loss_no_live'] = acc.get('loss_no_live', 0.)-float((aw*loss)[~has].sum())
            acc['loss_no_live_all_hidden'] = acc.get('loss_no_live_all_hidden', 0.)-float((aw*loss)[~has & hidden].sum())
            acc['loss_beyond_live_material'] = acc.get('loss_beyond_live_material', 0.)+float(
                (aw*np.maximum(-loss-ebar, 0.))[has].sum())
            acc['gain_no_live'] = acc.get('gain_no_live', 0.)+float((aw*gain)[~has].sum())
            ET0 = sim.erosion_totals.get('return_shortfall_thick_m3', 0.)
            out = orig(k, pl, wn, own, d, name)
            acc['shortfall_booked'] = acc.get('shortfall_booked', 0.)+sim.erosion_totals['return_shortfall_thick_m3']-ET0
            return out
        return orig(k, pl, wn, own, d, name)
    sim._tent_return = split
    for s in range(a.steps):
        sim.step()
        sim.state()
        if s % 20 == 19 or s == a.steps-1:
            print(sim.t, {k: f'{v/1e9:.4g} km3' for k, v in acc.items()}, flush=True)
    sb = acc['shortfall_booked']
    print('of the booked shortfall: no live contributor', acc['loss_no_live']/sb, '(all hidden',
          acc['loss_no_live_all_hidden']/sb, '); beyond the live material', acc['loss_beyond_live_material']/sb,
          '; the rest (the cap after rescaling)', 1-(acc['loss_no_live']+acc['loss_beyond_live_material'])/sb)


if __name__ == '__main__':
    main()
