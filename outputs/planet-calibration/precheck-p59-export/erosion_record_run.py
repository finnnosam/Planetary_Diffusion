"""Pass 59 pre-check: a plate run that records where the coarse erosion happens (output only; the model is unchanged).

Wraps scripts.simulate_plates. PlateSim._erode is wrapped to accumulate, over the last 40 Myr (t > -40), per world node:
the continental crust and sediment change of the erosion (e_cont, dsed: erosion < 0), this step's tectonic thickening
(self.tectonic, before erosion), and the surface the law saw (z_last); the last step's fields are kept too.
Written to <out>/seed<seed>/erosion_record.npz. Usage: same arguments as scripts.simulate_plates.
"""
import sys
from pathlib import Path

import numpy as np

from planet_tectonics import surface as surf

_orig = surf.SurfaceMixin._erode
REC = {}


def _erode(self, thick):
    e_cont, dsed = _orig(self, thick)
    if self.t > -40.+1e-9:
        a = REC.setdefault('acc', {k: np.zeros(self.g.size) for k in
                                   ('e_cont', 'dsed', 'tect_pos', 'tect', 'z0_sum', 'land_steps', 'orog_steps')})
        cont = self.world['cont']
        a['e_cont'] += e_cont
        a['dsed'] += dsed
        a['tect'] += self.tectonic
        a['tect_pos'] += np.clip(self.tectonic, 0, None)
        a['z0_sum'] += self.z_last
        a['land_steps'] += self.z_last > 0
        since = np.where(np.isnan(self.world['orog']), np.inf, self.t-self.world['orog'])
        a['orog_steps'] += cont & (since < 50.)
        REC['n'] = REC.get('n', 0)+1
        REC['last'] = {'e_cont': e_cont.copy(), 'dsed': dsed.copy(), 'tect': self.tectonic.copy(),
                       'z0': self.z_last.copy(), 'cont': cont.copy(), 'since': since.copy(),
                       'thick': thick.copy(), 'sed': self.world['sed'].copy()}
        if abs(self.t) < 1e-9:
            out = Path(REC['dir'])
            np.savez_compressed(out/'erosion_record.npz', steps=REC['n'], area=self.g.area,
                                **{f'acc_{k}': v for k, v in a.items()}, **{f'last_{k}': v for k, v in REC['last'].items()})
    return e_cont, dsed


surf.SurfaceMixin._erode = _erode

if __name__ == '__main__':
    import scripts.simulate_plates as sp
    args = sys.argv[1:]
    out = Path(args[args.index('--out')+1])
    seed = args[args.index('--seeds')+1]
    REC['dir'] = str(out/f'seed{seed}')
    (out/f'seed{seed}').mkdir(parents=True, exist_ok=True)
    sp.main(args)
