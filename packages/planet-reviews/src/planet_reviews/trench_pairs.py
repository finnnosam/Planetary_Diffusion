"""Ocean-ocean convergence plate pair by plate pair, planet against Earth (item B, T2; forty-fourth pass).

  python -m planet_reviews.trench_pairs outputs\\planet-calibration\\sweep-20261028-200644\\trench_polarity-node_trench_flux-band_arc_H_ocean_m-300_arc_accretion-False

Reads only saved snapshots and each seed's stats.json (event_list); changes nothing in the model. The object is
ocean-ocean convergence: on Earth Bird's (2003) SUB steps with oceanic crust 200 km to both sides (CRUST1.0 types, as
planet_reviews.trench_origin; 20.8k km) plus his OCB class (convergence without a slab, 17.4k km); on the planet the
contour segments typed subduction with oceanic crust on both sides (T2's intra-oceanic trench) plus slab-less
convergence (B-d3 runs only). Every segment gets:

- its closing rate (mm/yr);
- rho = |v_a - v_b| / max(|v_a|, |v_b|) at its midpoint, both plates' velocities in the no-net-rotation frame (rho
  near 0: the two plates move almost together; near 1 or above: independently or opposite);
- the floor age 100 km into each side (Earth: earth-seafloor-age-n217.npz, undated nodes left out; planet:
  ocean_age_myr), older and younger side;
- on the planet, the pair's lineage from the run's event list: the event kinds that first separated the two plates'
  lines of descent ('subduction initiation', 'rift', 'detach'), 'different start plates' (no common ancestor) and the time
  since they became separate plates; and since when the pair has had >= 500 km of intra-oceanic trench, read back
  every 10 Myr.

Planet distributions pool the snapshots at TIMES (length-weighted); pair tables are at 0 Myr. The same rho is printed
for every boundary class on both bodies, to tell a property of convergent pairs from one of all neighbours, and the
coherence of different plates' motions with distance: C(d) = sum v.v' / sum |v||v'| over node pairs on different plates
at great-circle distance d (both bodies NNR; 1 near-identical motion, 0 unrelated).
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_reviews import boundary_lengths as bl
from planet_reviews import force_balance_earth as fbe
from planet_reviews import trench_origin as to
from planet_tectonics import plates as P

TIMES = (0., -10., -20., -30.)
HISTORY = tuple(float(t) for t in range(0, -251, -10))
RHO_CO = .5            # "moving almost together"
YOUNG_MYR = 20.
SIDE_KM = 100.
R_KM = 6371.
SINCE_BINS = ((0, 10), (10, 30), (30, 100), (100, 1e9))
COH_BINS_KM = (0, 500, 1000, 2000, 3000, 5000, 8000)


def coherence(g, pid, vel, rng, n_ref=1500, n_cmp=40000):
    """C(d) over node pairs on different plates, area-weighted samples."""
    p = g.area/g.area.sum()
    i = rng.choice(g.size, n_ref, p=p); j = rng.choice(g.size, n_cmp, p=p)
    num = np.zeros(len(COH_BINS_KM)-1); den = np.zeros_like(num)
    vj, xj, pj = vel[j], g.xyz[j], pid[j]
    nj = np.linalg.norm(vj, axis=1)
    edges = np.asarray(COH_BINS_KM, float)/R_KM
    for a in i:
        d = np.arccos(np.clip(xj@g.xyz[a], -1, 1))
        ok = pj != pid[a]
        b = np.searchsorted(edges, d[ok])-1
        k = (b >= 0) & (b < num.size)
        num += np.bincount(b[k], weights=(vj[ok][k]@vel[a]), minlength=num.size)
        den += np.bincount(b[k], weights=nj[ok][k]*np.linalg.norm(vel[a]), minlength=num.size)
    return (num/np.maximum(den, 1e-30)).tolist()


def coh_txt(c):
    return ', '.join(f'{lo}-{hi} km {v:.2f}' for lo, hi, v in zip(COH_BINS_KM[:-1], COH_BINS_KM[1:], c))


def wq(v, w, qs=(.1, .5, .9)):
    v, w = np.asarray(v, float), np.asarray(w, float)
    k = np.isfinite(v)
    return to.wpct(v[k], w[k], qs)


def share(mask, w):
    w = np.asarray(w, float)
    return float(w[mask].sum()/w.sum()) if w.sum() > 0 else float('nan')


def side_points(mid, nrm, km):
    pa = mid+km/R_KM*nrm; pb = mid-km/R_KM*nrm
    return pa/np.linalg.norm(pa, axis=1, keepdims=True), pb/np.linalg.norm(pb, axis=1, keepdims=True)


def rho(va, vb):
    return np.linalg.norm(va-vb, axis=1)/np.maximum(np.maximum(np.linalg.norm(va, axis=1), np.linalg.norm(vb, axis=1)), 1e-9)


# ----------------------------------------------------------------------------------------------------------- Earth

def earth_segments():
    """Bird's steps with class, length (km), closing rate (mm/yr), rho, both side ages and pair names."""
    g, pid, names, _ = fbe.earth_plates(217)
    w = fbe.earth_omega(g, pid, names)                              # rad/s, NNR
    st = fbe.earth_steps(names, fbe._rings())
    x, nA = st['x'], st['nA']
    ms = 1e3*fbe.YR                                                 # m/s -> mm/yr
    va = np.cross(np.stack([w[int(k)] for k in st['A']]), x)*fbe.R*ms
    vb = np.cross(np.stack([w[int(k)] for k in st['B']]), x)*fbe.R*ms
    types = np.array([l.split() for l in (DATA/'crust1'/'CNtype1-1.txt').read_text().splitlines() if l.strip()])

    def oceanic(dist_km):                                           # as trench_origin.earth_reference
        ok = np.ones(len(x), bool)
        for sg in (1, -1):
            for d in (dist_km/2, dist_km):
                q = x+sg*d/R_KM*nA; q /= np.linalg.norm(q, axis=1, keepdims=True)
                lat = np.degrees(np.arcsin(np.clip(q[:, 2], -1, 1))); lon = np.degrees(np.arctan2(q[:, 1], q[:, 0]))
                r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
                ok &= np.isin(types[r, c], to.OCEANIC_TYPES)
        return ok
    za = np.load(OUTPUTS/'planet-reference'/'earth-seafloor-age-n217.npz')
    age = np.where(za['dated_share'] > .5, za['age_myr'], np.nan)
    pa, pb = side_points(x, nA, SIDE_KM)                            # pa on B's side (nA points out of A)
    ageB, ageA = age[fields.nearest_nodes(g, pa)], age[fields.nearest_nodes(g, pb)]
    cls = np.asarray(st['cls'])
    io = (cls == 'SUB') & oceanic(200)
    nm = np.array(names)
    wn = np.stack([w[int(k)] for k in pid])
    coh = coherence(g, pid, np.cross(wn, g.xyz), np.random.default_rng(0))
    return {'coh': coh, 'L': st['L']/1e3, 'cls': cls, 'close': -np.asarray(st['div'], float), 'rho': rho(va, vb),
            'dv': np.linalg.norm(va-vb, axis=1), 'old': np.fmax(ageA, ageB), 'young': np.fmin(ageA, ageB),
            'io': io, 'ocb': cls == 'OCB', 'pair': np.array([f'{a}-{b}' if a < b else f'{b}-{a}'
                                                             for a, b in zip(nm[st['A']], nm[st['B']])])}


# ---------------------------------------------------------------------------------------------------------- planet

def lineage(events):
    """parent[id] = (parent id, kind, t) for plates created by split events."""
    par = {}
    for e in events:
        if e['kind'] in ('rift', 'subduction initiation', 'detach') and 'new_plate' in e:
            par[int(e['new_plate'])] = (int(e['plate']), e['kind'], float(e['t']))
    return par


def pair_relation(a, b, par, t0=-250.):
    """(kinds that first separated the two lines, time since they have been separate plates (Myr before the
    snapshot is applied by the caller: returns the separation time), 'different start plates' when the roots differ)."""
    def chain(k):
        out = [(k, None, None)]
        while k in par:
            p, kind, t = par[k]
            out[-1] = (out[-1][0], kind, t)
            out.append((p, None, None)); k = p
        return out                                                  # [(id, kind that made it, t), ..., (root, None, None)]
    ca, cb = chain(a), chain(b)
    ia = {k: i for i, (k, _, _) in enumerate(ca)}
    for j, (k, _, _) in enumerate(cb):
        if k in ia:
            i = ia[k]
            first = [c for c in (ca[i-1] if i > 0 else None, cb[j-1] if j > 0 else None) if c is not None]
            return '+'.join(sorted({c[1] for c in first})), min(c[2] for c in first)
    return 'different start plates', t0


def planet_snapshot(g, f, par):
    z = np.load(str(f)+'.npz'); plates = json.loads(Path(str(f)+'-plates.json').read_text())
    t = float(z['time_myr'])
    pid = z['plate_id'].astype(np.int64); cont = z['continental'].astype(bool)
    s = bl.typed_segments(g, pid, plates, cont)
    L, typ, cv, na, nb, ka, kb, mid = s['L'], s['typ'], s['conv'], s['na'], s['nb'], s['ka'], s['kb'], s['mid']
    oo = ~(cont[na] | cont[nb])
    bts = z['boundary_type']
    side_sub = (bts[na] == P.SUBDUCTION) | (bts[nb] == P.SUBDUCTION)
    noslab = (typ == 2) & ((bts[na] == P.CONVERGENT) | (bts[nb] == P.CONVERGENT)) & ~side_sub
    io = (typ == 2) & oo & ~noslab
    ocb = noslab & oo
    om = {int(q['id']): np.asarray(q['axis'], float)*q['speed_cm_yr']*1e-2/g.radius for q in plates}   # rad/yr
    w = fbe.remove_net_rotation(g, pid, om)
    W = lambda ks: np.stack([w.get(int(k), np.zeros(3)) for k in ks])
    va = np.cross(W(ka), mid)*g.radius*1e3; vb = np.cross(W(kb), mid)*g.radius*1e3      # mm/yr
    d = g.xyz[na]-g.xyz[nb]; d -= (d*mid).sum(1, keepdims=True)*mid
    d /= np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-15)
    pa, pb = side_points(mid, d, SIDE_KM)
    age = np.where(cont, np.nan, z['ocean_age_myr'].astype(float))
    A, B = age[fields.nearest_nodes(g, pa)], age[fields.nearest_nodes(g, pb)]
    keys = np.array([f'{min(a, b)}-{max(a, b)}' for a, b in zip(ka, kb)])
    rel = {}
    for k in np.unique(keys[io | ocb]):
        a, b = (int(x) for x in k.split('-'))
        kind, ts = pair_relation(a, b, par)
        rel[k] = (kind, t-ts)
    kinds = np.array([rel.get(k, ('', np.nan))[0] for k in keys])
    since = np.array([rel.get(k, ('', np.nan))[1] for k in keys])
    cls = np.select([typ == 1, typ == 4, io, ocb, typ == 3, typ == 2], ['ridge', 'transform', 'io', 'ocb', 'collision',
                                                                     'cont-sub'], 'none')
    wn = np.stack([w.get(int(k), np.zeros(3)) for k in pid])
    coh = coherence(g, pid, np.cross(wn, g.xyz), np.random.default_rng(0))
    return {'coh': coh, 't': t, 'L': L, 'cls': cls, 'close': cv, 'rho': rho(va, vb), 'dv': np.linalg.norm(va-vb, axis=1),
            'old': np.fmax(A, B), 'young': np.fmin(A, B), 'pair': keys, 'kind': kinds, 'since': since,
            'io': io, 'ocb': ocb, 'plates': {int(q['id']): q for q in plates}}


def io_history(g, sd):
    """Intra-oceanic trench length by pair every 10 Myr (k: pair -> {t: km})."""
    h = defaultdict(dict)
    for t in HISTORY:
        f = to.snap_path(sd, t)
        if not Path(str(f)+'.npz').exists():
            continue
        z = np.load(str(f)+'.npz'); plates = json.loads(Path(str(f)+'-plates.json').read_text())
        pid = z['plate_id'].astype(np.int64); cont = z['continental'].astype(bool)
        s = bl.typed_segments(g, pid, plates, cont)
        bts = z['boundary_type']
        na, nb = s['na'], s['nb']
        noslab = (s['typ'] == 2) & ((bts[na] == P.CONVERGENT) | (bts[nb] == P.CONVERGENT)) & \
            ~((bts[na] == P.SUBDUCTION) | (bts[nb] == P.SUBDUCTION))
        io = (s['typ'] == 2) & ~(cont[na] | cont[nb]) & ~noslab
        for a, b, l in zip(s['ka'][io], s['kb'][io], s['L'][io]):
            k = f'{min(a, b)}-{max(a, b)}'
            h[k][t] = h[k].get(t, 0.)+float(l)
    return h


def io_since(hist, key, km=500.):
    """Earliest time of the unbroken run of 10 Myr snapshots, back from 0, with >= km of the pair's trench."""
    since = None
    for t in HISTORY:
        if hist.get(key, {}).get(t, 0.) >= km:
            since = t
        else:
            break
    return since


# ---------------------------------------------------------------------------------------------------------- report

def summary(seg, m):
    L = seg['L'][m]
    return {'km': float(L.sum()), 'close_p10_p50_p90': wq(seg['close'][m], L), 'slow20': share(seg['close'][m] < 20, L),
            'rho_p10_p50_p90': wq(seg['rho'][m], L), 'rho_lt': share(seg['rho'][m] < RHO_CO, L),
            'dv_p50': wq(seg['dv'][m], L, (.5,))[0], 'old_p10_p50_p90': wq(seg['old'][m], L),
            'old_young': share(seg['old'][m] < YOUNG_MYR, L*np.isfinite(seg['old'][m])),
            'young_p50': wq(seg['young'][m], L, (.5,))[0], 'age_dated': share(np.isfinite(seg['old'][m]), L)}


def fmt(name, q, per=1.):
    return (f"  {name:<34s} {q['km']/per/1e3:6.1f}k  close {'/'.join(f'{x:.0f}' for x in q['close_p10_p50_p90'])}"
            f" (<20: {q['slow20']:.2f})  rho {'/'.join(f'{x:.2f}' for x in q['rho_p10_p50_p90'])}"
            f" (<{RHO_CO}: {q['rho_lt']:.2f})  |dv| p50 {q['dv_p50']:.0f}  older side age"
            f" {'/'.join(f'{x:.0f}' for x in q['old_p10_p50_p90'])} (<{YOUNG_MYR:.0f}: {q['old_young']:.2f}),"
            f" younger p50 {q['young_p50']:.0f}, dated {q['age_dated']:.2f}")


def cat(segs):
    out = {}
    for k in ('L', 'cls', 'close', 'rho', 'dv', 'old', 'young', 'io', 'ocb'):
        out[k] = np.concatenate([s[k] for s in segs])
    for k in ('kind', 'since'):
        if k in segs[0]:
            out[k] = np.concatenate([s[k] for s in segs])
    return out


def seed_dirs(dirs):
    return to.seed_dirs(dirs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='+', help='combination folders (all seeds) or seed folders')
    ap.add_argument('--no-history', action='store_true', help='skip the 10 Myr read-back (faster)')
    ap.add_argument('--json', help='write the numbers here')
    a = ap.parse_args(argv)
    out = {}
    E = earth_segments()
    print('EARTH (Bird 2003 steps, NNR poles; CRUST1.0 types; seafloor age 100 km to each side)')
    eoo = E['io'] | E['ocb']
    out['earth'] = {'oo_convergence': summary(E, eoo), 'sub_io': summary(E, E['io']), 'ocb': summary(E, E['ocb'])}
    print(fmt('ocean-ocean convergence (SUB io+OCB)', out['earth']['oo_convergence']))
    print(fmt('  SUB intra-oceanic', out['earth']['sub_io']))
    print(fmt('  OCB (no slab)', out['earth']['ocb']))
    print('  rho by class, all steps (length-weighted p50; share < %.1f):' % RHO_CO)
    out['earth']['rho_by_class'] = {}
    for c in ('OSR', 'CRB', 'OTF', 'CTF', 'SUB', 'OCB', 'CCB'):
        m = E['cls'] == c
        out['earth']['rho_by_class'][c] = (wq(E['rho'][m], E['L'][m], (.5,))[0], share(E['rho'][m] < RHO_CO, E['L'][m]),
                                           float(E['L'][m].sum()), wq(E['dv'][m], E['L'][m], (.5,))[0])
    print('   ' + '; '.join(f'{c} {v[0]:.2f} ({v[1]:.2f}, {v[2]/1e3:.1f}k, |dv| {v[3]:.0f})' for c, v in out['earth']['rho_by_class'].items()))
    out['earth']['coherence'] = E['coh']
    print('  coherence of different plates\' motions by distance: ' + coh_txt(E['coh']))
    print('  pairs (k km, closing p50 mm/yr, rho p50, older-side age p50 Myr):')
    for lab, m in (('SUB io', E['io']), ('OCB', E['ocb'])):
        pairs = defaultdict(float)
        for p, l in zip(E['pair'][m], E['L'][m]):
            pairs[p] += l
        top = sorted(pairs.items(), key=lambda x: -x[1])[:10]
        txt = []
        for p, l in top:
            k = m & (E['pair'] == p)
            txt.append(f"{p} {l/1e3:.1f} ({wq(E['close'][k], E['L'][k], (.5,))[0]:.0f}, {wq(E['rho'][k], E['L'][k], (.5,))[0]:.2f},"
                       f" {wq(E['old'][k], E['L'][k], (.5,))[0]:.0f})")
        print(f'   {lab}: ' + '; '.join(txt))

    g = CubeSphere(217)
    allseg, rows0 = [], []
    for sd in seed_dirs(a.dirs):
        st = json.loads((Path(sd)/'stats.json').read_text())
        par = lineage(st.get('event_list', []))
        segs = []
        for t in TIMES:
            f = to.snap_path(sd, t)
            if Path(str(f)+'.npz').exists():
                segs.append(planet_snapshot(g, f, par))
        allseg += segs
        hist = {} if a.no_history else io_history(g, sd)
        s0 = segs[0]
        pairs = defaultdict(float)
        for p, l in zip(s0['pair'][s0['io']], s0['L'][s0['io']]):
            pairs[p] += l
        pr = []
        for p, l in sorted(pairs.items(), key=lambda x: -x[1]):
            k = s0['io'] & (s0['pair'] == p)
            x, y = (int(v) for v in p.split('-'))
            cf = [s0['plates'][i]['continental_km2']/s0['plates'][i]['area_km2'] for i in (x, y)]
            ar = [s0['plates'][i]['area_km2']/1e6 for i in (x, y)]
            pr.append({'pair': p, 'km': l, 'close_p50': wq(s0['close'][k], s0['L'][k], (.5,))[0],
                       'rho_p50': wq(s0['rho'][k], s0['L'][k], (.5,))[0], 'old_p50': wq(s0['old'][k], s0['L'][k], (.5,))[0],
                       'kind': str(s0['kind'][k][0]), 'separate_myr': float(s0['since'][k][0]),
                       'io_since': io_since(hist, p) if hist else None, 'area_mkm2': ar, 'cf': cf})
        rows0.append((str(sd), sum(pairs.values()), pr))
    P_ = cat(allseg)
    nsnap = len(allseg)
    poo = P_['io'] | P_['ocb']
    print(f'\nPLANET ({len(rows0)} seeds x {nsnap // max(len(rows0), 1)} snapshots {TIMES}; km per snapshot)')
    out['planet'] = {'oo_convergence': summary(P_, poo), 'io': summary(P_, P_['io']), 'ocb': summary(P_, P_['ocb']),
                     'snapshots': nsnap}
    print(fmt('ocean-ocean convergence (io+ocb)', out['planet']['oo_convergence'], nsnap))
    print(fmt('  io (T2)', out['planet']['io'], nsnap))
    if P_['ocb'].any():
        print(fmt('  ocb (no slab)', out['planet']['ocb'], nsnap))
    print('  rho by class (length-weighted p50; share < %.1f; km per snapshot):' % RHO_CO)
    out['planet']['rho_by_class'] = {}
    for c in ('ridge', 'transform', 'io', 'ocb', 'cont-sub', 'collision'):
        m = P_['cls'] == c
        if m.any():
            out['planet']['rho_by_class'][c] = (wq(P_['rho'][m], P_['L'][m], (.5,))[0], share(P_['rho'][m] < RHO_CO, P_['L'][m]),
                                                float(P_['L'][m].sum()/nsnap), wq(P_['dv'][m], P_['L'][m], (.5,))[0])
    print('   ' + '; '.join(f'{c} {v[0]:.2f} ({v[1]:.2f}, {v[2]/1e3:.1f}k, |dv| {v[3]:.0f})' for c, v in out['planet']['rho_by_class'].items()))
    coh = np.mean([s['coh'] for s in allseg], axis=0).tolist()
    out['planet']['coherence'] = coh
    print('  coherence of different plates\' motions by distance (mean of snapshots): ' + coh_txt(coh))
    # lineage
    print('  io trench by lineage of the pair (share of io length; closing p50; rho p50; older-side age p50):')
    m = P_['io']
    out['planet']['lineage'] = {}
    for kind in sorted(set(P_['kind'][m])):
        k = m & (P_['kind'] == kind)
        q = {'share': share(k[m], P_['L'][m]), 'close_p50': wq(P_['close'][k], P_['L'][k], (.5,))[0],
             'rho_p50': wq(P_['rho'][k], P_['L'][k], (.5,))[0], 'old_p50': wq(P_['old'][k], P_['L'][k], (.5,))[0],
             'since': {f'{lo}-{hi if hi < 1e8 else "start"}': share(k[m] & (P_['since'][m] >= lo) & (P_['since'][m] < hi), P_['L'][m])
                       for lo, hi in SINCE_BINS}}
        out['planet']['lineage'][kind] = q
        print(f"   {kind:<34s} {q['share']:.2f}  close {q['close_p50']:.0f}  rho {q['rho_p50']:.2f}  old {q['old_p50']:.0f};"
              f" separate for (Myr) " + ', '.join(f'{b}: {v:.2f}' for b, v in q['since'].items()))
    print('\n  pairs at 0 Myr (k km; closing p50; rho p50; older-side age p50; lineage, separate for; io since; areas Mkm2;'
          ' continental shares):')
    out['planet']['pairs0'] = []
    for sd, tot, pr in rows0:
        print(f'   {Path(sd).name}: io {tot/1e3:.1f}k')
        for q in pr[:8]:
            print(f"     {q['pair']:>7s} {q['km']/1e3:4.1f}  {q['close_p50']:3.0f}  {q['rho_p50']:.2f}  {q['old_p50']:4.0f}  "
                  f"{q['kind']}, {q['separate_myr']:.0f}; since {q['io_since'] if q['io_since'] is not None else '-'};"
                  f" {q['area_mkm2'][0]:.0f}/{q['area_mkm2'][1]:.0f}; {q['cf'][0]:.2f}/{q['cf'][1]:.2f}")
        out['planet']['pairs0'].append({'seed': sd, 'io_km': tot, 'pairs': pr})
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1, default=float))
    return out


if __name__ == '__main__':
    main()
