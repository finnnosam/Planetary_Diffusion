"""Rifted-margin production through a run against Earth's since 250 Ma (sixty-fifth pass, journal 87; candidate K).

  python -m planet_reviews.rift_production <sweep or combination or seed folders...>

At each snapshot every 10 Myr: the continental nodes beside ocean floor <= 10 Myr old, leaving out those within
300 km of a subduction boundary (their young floor is a back-arc or a young plate at a trench, not a breakup). Floor
that young beside a margin was made in the last 10 Myr, so length / 10 Myr is the production of rifted margin over
that interval, each margin counted at its birth (no survival bias); both conjugate margins count. Length = nodes x
graph spacing. Continent: the continental flag (the counterpart of the age grid's non-ocean area); the label-free
reading (crust + sediment >= 20 km) is printed beside it.

Earth (planet_tectonics.breakup.EARTH_RIFT_PRODUCTION, the table split_rule 'earth' uses; the same instrument on
the same n = 217 graph, pass 65): Mueller et al. (2019) paleo-age grids
every 10 Myr, subduction zones resolved from the same model with pygplates (`earth_rift_production.py` in
outputs/planet-calibration/precheck-p65-towers/). Mean 1,553 km/Myr over 250 Myr (7.4 km/Myr per Mkm^2 of non-ocean
area), 388k km made in all: low at 250-170 Ma (Pangea whole), the breakup peak at 160-110 Ma (2,300-4,300 km/Myr),
1,000-2,900 km/Myr since.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_tectonics.breakup import EARTH_RIFT_PRODUCTION as EARTH   # (t Ma, km/Myr over (t, t+10] Ma, Mkm^2)
EARTH_MEAN_KM_MYR = float(np.mean([r[1] for r in EARTH]))
EARTH_MEAN_PER_MKM2 = float(np.mean([r[1]/r[2] for r in EARTH]))
SUBDUCTION = 2
NEAR_KM = 300.
YOUNG_MYR = 10.


def production(z, g, params, sp):
    """(rifted km, all young-edge km, continental area Mkm^2, label-free rifted km) for one snapshot."""
    from planet_landscape import graph_landscape as gl
    from planet_reviews.crust_review import column_thickness_m, label_free_continent
    nb = g.nbr; ok = nb >= 0
    nbv = lambda m: np.any(np.where(ok, m[np.where(ok, nb, 0)], False), 1)
    cont = z['continental'].astype(bool)
    age = z['ocean_age_myr'].astype(np.float64)
    young = (~cont) & np.isfinite(age) & (age <= YOUNG_MYR)
    sz = z['boundary_type'] == SUBDUCTION
    near = gl.distance_to(sz, nb, g.dist, (NEAR_KM+100)*1e3) <= NEAR_KM*1e3 if sz.any() else np.zeros(g.size, bool)
    ye = cont & nbv(young)
    lf = label_free_continent(column_thickness_m(z, params))
    ye_lf = lf & nbv(young & ~lf)
    return ((ye & ~near).sum()*sp/1e3, ye.sum()*sp/1e3, g.area[cont].sum()/1e12, (ye_lf & ~near).sum()*sp/1e3)


def read_seed(sd, g, sp):
    from planet_reviews.crust_review import run_params
    p = run_params(sd)
    out = {}
    for f in sorted(Path(sd, 'snapshots').glob('t*.npz')):
        z = np.load(f)
        t = float(z['time_myr'])
        if abs(t/10-round(t/10)) > 1e-6:
            continue
        r, a, A, lf = production(z, g, p, sp)
        out[round(t)] = {'rifted_km': r, 'young_edge_km': a, 'cont_mkm2': A, 'rifted_lf_km': lf}
    st = json.loads(Path(sd, 'stats.json').read_text()) if Path(sd, 'stats.json').exists() else {}
    ev = st.get('events', {}) if isinstance(st.get('events'), dict) else {}
    return out, ev


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('folders', nargs='+')
    ap.add_argument('--json')
    a = ap.parse_args(argv)
    from planet_grid.sphere_graph import CubeSphere
    g = CubeSphere(217); sp = float(np.sqrt(g.area.mean()))
    print(f"Earth (Mueller et al. 2019, same instrument): mean {EARTH_MEAN_KM_MYR:.0f} km/Myr "
          f"({EARTH_MEAN_PER_MKM2:.2f} per Mkm^2 of continent); 250-150 Ma {np.mean([r[1] for r in EARTH[15:]]):.0f}, "
          f"150-50 Ma {np.mean([r[1] for r in EARTH[5:15]]):.0f}, 50-0 Ma {np.mean([r[1] for r in EARTH[:5]]):.0f}")
    seeds = []
    for f in map(Path, a.folders):
        seeds += [f] if (f/'snapshots').is_dir() else sorted(p for p in f.rglob('seed*') if (p/'snapshots').is_dir())
    combos = {}
    for s in seeds:
        combos.setdefault(s.parent, []).append(s)
    res = {}
    for c, ss in combos.items():
        per = {s.name: read_seed(s, g, sp)[0] for s in ss}
        ts = sorted(set.intersection(*[set(v) for v in per.values()]))
        ts = [t for t in ts if t > -250]                 # the start snapshot holds no 10 Myr of history
        P = np.array([[per[s][t]['rifted_km']/YOUNG_MYR for t in ts] for s in per])
        Pa = np.array([[per[s][t]['rifted_km']/YOUNG_MYR/per[s][t]['cont_mkm2'] for t in ts] for s in per])
        PL = np.array([[per[s][t]['rifted_lf_km']/YOUNG_MYR for t in ts] for s in per])
        Y = np.array([[per[s][t]['young_edge_km']/YOUNG_MYR for t in ts] for s in per])
        tt = np.array(ts)
        early, mid, late = tt <= -150, (tt > -150) & (tt <= -50), tt > -50
        print(f"{c.name} ({len(per)} seeds): rifted production mean {P.mean():.0f} km/Myr "
              f"({Pa.mean():.2f} per Mkm^2; Earth {EARTH_MEAN_PER_MKM2:.2f}), x{Pa.mean()/EARTH_MEAN_PER_MKM2:.2f} of Earth per area; "
              f"by seed {' / '.join(f'{v:.0f}' for v in P.mean(1))}; -250..-150 {P[:, early].mean():.0f}, -150..-50 "
              f"{P[:, mid].mean():.0f}, -50..0 {P[:, late].mean():.0f}; label-free {PL.mean():.0f}; all young edge {Y.mean():.0f}")
        print('    t (Myr): ' + ' '.join(f'{t:5d}' for t in ts))
        print('    planet : ' + ' '.join(f'{v:5.0f}' for v in P.mean(0)))
        er = {-r[0]: r[1] for r in EARTH}
        print('    Earth  : ' + ' '.join(f"{er.get(t, np.nan):5.0f}" for t in ts))
        res[c.name] = {'t': ts, 'prod_km_myr': P.mean(0).tolist(), 'per_mkm2_mean': float(Pa.mean()),
                       'mean_km_myr': float(P.mean()), 'by_seed': P.mean(1).tolist()}
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
