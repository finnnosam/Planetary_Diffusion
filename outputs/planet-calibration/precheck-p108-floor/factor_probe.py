"""Pass 109 (journal 147): the remap's rescaling factors. From the warm state, K steps under a rule; every
`_conserved` call's gain and loss factors (world volume / plate volume) and the plates' thickness range after each step.

  python factor_probe.py <warm dir> <rule> [K]
"""
import copyreg
import importlib
import pickle
import sys
import types

import numpy as np

copyreg.pickle(types.ModuleType, lambda m: (importlib.import_module, (m.__name__,)))
from planet_tectonics.params import PARAMS  # noqa: E402

warm, rule = sys.argv[1], sys.argv[2]
K = int(sys.argv[3]) if len(sys.argv) > 3 else 5
sim = pickle.load(open(f'{warm}/warm.pkl', 'rb'))
for k, v in PARAMS.items():
    sim.p.setdefault(k, v)
sim.p['remap_rule'] = rule
sim.state()
rec = []
orig = sim._conserved


def conserved(dj, aj, dw, aw, name):
    dj64, dw64 = np.asarray(dj, np.float64), np.asarray(dw, np.float64)
    row = {'name': name, 'n': int(dj64.size)}
    for sg, lab in ((1., 'gain'), (-1., 'loss')):
        Pj = float((aj*dj64)[sg*dj64 > 0].sum())
        Pw = float((aw*dw64)[sg*dw64 > 0].sum())
        row[lab] = Pw/Pj if Pj != 0 else (np.inf if Pw != 0 else 1.)
        row[lab+'_Pw_km3'] = Pw/1e9
    row['net_ratio'] = float((aw*dw64).sum())/float((aj*dj64).sum()) if float((aj*dj64).sum()) != 0 else np.nan
    rec.append(row)
    return orig(dj, aj, dw, aw, name)


sim._conserved = conserved
for i in range(K):
    rec.clear()
    sim.step()
    sim.state()
    th = np.concatenate([pl.thick[pl.mask].astype(np.float64) for pl in sim.plates])
    sd = np.concatenate([pl.sed[pl.mask].astype(np.float64) for pl in sim.plates])
    big = [r for r in rec if max(abs(r['gain']), abs(r['loss'])) > 3 or min(r['gain'], r['loss']) < 0.33]
    g = np.array([r['gain'] for r in rec]); l = np.array([r['loss'] for r in rec])
    print(f"{rule} t={sim.t:+.0f}: {len(rec)} calls; gain factor median {np.median(g):.3f} max {np.max(g):.3g}; "
          f"loss factor median {np.median(l):.3f} max {np.max(l):.3g}; calls outside 1/3..3: {len(big)}; "
          f"plate thick min {th.min():.0f} max {th.max():.0f} m; sed max {sd.max():.0f} m", flush=True)
    for r in sorted(big, key=lambda r: -max(r['gain'], r['loss']))[:4]:
        print(f"   {r['name']} n {r['n']}: gain x{r['gain']:.3g} (world {r['gain_Pw_km3']:.3g} km3), loss x{r['loss']:.3g} "
              f"(world {r['loss_Pw_km3']:.3g} km3), net ratio {r['net_ratio']:.3g}", flush=True)
