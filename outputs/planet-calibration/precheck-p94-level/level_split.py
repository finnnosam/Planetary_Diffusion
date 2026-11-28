"""Pass 94 (diagnosis): the continents' level on the datum, split by crust bin and term, planet against Earth with
alike statistics on the same n = 217 nodes. Run from the repo root.

Domain: label-free continent, total crust + bulk sediment >= 20 km (planet: the continental flag as well; Earth:
CRUST1.0 non-oceanic type, ice boxes out). Bins of total thickness 20-25 ... 45-55, 55-80 km.
Per bin, area-weighted: share of the surface, share under the datum, z mean and p25 / p50 / p75 on the datum, and the
mean of each term of the model's own column law (common.py): Airy on the crystalline crust, densification, sediment
column, lithosphere C, dynamic topography, rest. Means add, so planet - Earth splits exactly into the terms' gaps.
Earth: C from Hoggard 2020 L (c_seis, clipped to the frame as the model clips) and from the class age (c_age) beside it;
dyn from Holdt et al. 2022 degree 1-30 (the continental convention unconfirmed, journal 73), the rest without it.
Planet: today's control (sweep-20261119-182758 arm C), 0 Myr, six seeds (sd across seeds).

RULE fixed before the run (pass 94 orders pass 95's candidates):
  weight w_b = Earth's area share of bin b among 30-45 km (where the drowned continent lies: 30-35 / 35-40 / 40-45);
  gap_t = sum_b w_b (planet - Earth)(mean term t), for t in C, sediment (sedc + the sediment's part of dens is not
  separable: dens read with it), dyn. The candidates rank by -gap_t (the term holding the planet lowest first):
  (a) the L-T coupling <- C; (b) the drowned-continent sediment column <- sediment + dens; (c) dyn's placement <- dyn.
  A term whose gap is within its seed sd/sqrt(6)*2 of zero does not rank. Airy and rest are read, not candidates
  (Airy: the crust's distribution, passed; rest: flexure / water load, a remainder)."""
import sys, numpy as np
sys.path.insert(0, 'outputs/planet-calibration/precheck-p82-stepback')
import common as C
ARM = sys.argv[1] if len(sys.argv) > 1 else 'outputs/planet-calibration/sweep-20261119-182758/backarc_extension-True_rollback_rule-fixed_rollback_mm_yr-9'
B = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (55, 80)]
TERMS = ['airy', 'dens', 'sedc', 'C', 'dyn', 'rest']
e = C.earth(C.run_params(f'{ARM}/seed{C.SEEDS[0]}'))
from planet_tectonics import lithosphere as lt
try:
    import h5py
    D = 'data/residual_topography/Residual depthtopographic measurements/Residual depthtopographic measurements/Holdt_et_al_2022_SupportingInformation/'
    with h5py.File(D+'degree_1_to_30.grd', 'r') as f:
        lon, lat, v = f['lon'][:], f['lat'][:], f['z'][:]
    i = np.clip(np.round((e['lat']-lat[0])/(lat[1]-lat[0])).astype(int), 0, lat.size-1)
    j = np.round(((e['lon'] % 360)-lon[0])/(lon[1]-lon[0])).astype(int) % lon.size
    edyn = np.nan_to_num(v[i, j]*1e3)
except Exception as x:
    print('Holdt field not read:', x); edyn = np.zeros(e['z'].size)
ED = dict(airy=e['airy'], dens=e['dens'], sedc=e['sedc'], C=e['c_seis_clip'], dyn=edyn)
ED['rest'] = e['z']-sum(ED[t] for t in TERMS[:-1])
ECA = e['c_age']
dom_e = e['cont'] & ~e['ice'] & np.isfinite(e['L'])
ae = e['a']/e['a'].sum()

def stats(z, a, m, terms, sed, L):
    a_m = a[m]; out = dict(share=a_m.sum(), under=a_m[z[m] < 0].sum()/max(a_m.sum(), 1e-30))
    out['mean'] = C.wmean(z[m], a_m); out.update(zip(('p25', 'p50', 'p75'), C.wq(z[m], a_m, [.25, .5, .75])))
    for t in TERMS: out[t] = C.wmean(terms[t][m], a_m)
    out['sed'] = C.wmean(sed[m], a_m); out['L'] = C.wmean(L[m], a_m)
    return out

E = {}
for lo, hi in B:
    m = dom_e & (e['tot'] >= lo) & (e['tot'] < hi)
    E[(lo, hi)] = stats(e['z'], ae, m, ED, e['S'], e['L']); E[(lo, hi)]['C_age'] = C.wmean(ECA[m], ae[m])
Pl = {b: [] for b in B}; land = []
for s in C.SEEDS:
    P = C.planet(s, 't+0000.0', ARM); a = P['a']/P['a'].sum()
    T = dict(airy=P['airy'], dens=P['dens'], sedc=P['sedc'], C=P['cterm'], dyn=P['dyn'], rest=P['rest'])
    for lo, hi in B:
        m = P['cont'] & (P['tot'] >= lo) & (P['tot'] < hi)
        Pl[(lo, hi)].append(stats(P['zd'], a, m, T, P['sh'], np.nan_to_num(P['L'])))
keys = ['share', 'under', 'mean', 'p25', 'p50', 'p75'] + TERMS + ['sed', 'L']
PM = {b: {k: np.mean([r[k] for r in Pl[b]]) for k in keys} for b in B}
PS = {b: {k: np.std([r[k] for r in Pl[b]], ddof=1) for k in keys} for b in B}
print(ARM)
print('\n1. Per bin (total crust + bulk sediment): surface share | share under the datum | z mean / p25 / p50 / p75 on the datum (m) | mean sediment (m) | mean L (km)')
for b in B:
    p, q = PM[b], E[b]
    print(f'  {b[0]:2d}-{b[1]:2d} km  planet {p["share"]:.3f} | {p["under"]:.2f} | {p["mean"]:6.0f} ({PS[b]["mean"]:3.0f}) / {p["p25"]:6.0f} / {p["p50"]:6.0f} / {p["p75"]:6.0f} | {p["sed"]:5.0f} | {p["L"]:4.0f}')
    print(f'            Earth  {q["share"]:.3f} | {q["under"]:.2f} | {q["mean"]:6.0f}       / {q["p25"]:6.0f} / {q["p50"]:6.0f} / {q["p75"]:6.0f} | {q["sed"]:5.0f} | {q["L"]:4.0f}')
print('\n2. Mean of each term (m): planet / Earth / gap (planet - Earth; seed sd of the planet term). Earth C from Hoggard L (C by class age in brackets); Earth dyn Holdt 1-30')
for b in B:
    p, q, s = PM[b], E[b], PS[b]
    print(f'  {b[0]:2d}-{b[1]:2d} km  z gap {p["mean"]-q["mean"]:+5.0f} = ' + ' + '.join(f'{t} {p[t]-q[t]:+5.0f}' for t in TERMS))
    print('            ' + ' | '.join(f'{t} {p[t]:6.0f} / {q[t]:6.0f} ({s[t]:3.0f})' for t in TERMS) + f' | C_age {q["C_age"]:5.0f}')
W = np.array([E[b]['share'] for b in B[2:5]]); W = W/W.sum()
print(f'\n3. RULE: Earth-share weights over 30-35 / 35-40 / 40-45 km: ' + ' / '.join(f'{w:.2f}' for w in W))
gap = {}; sem = {}
for t in TERMS + ['mean']:
    seeds = np.array([[r[t] for r in Pl[b]] for b in B[2:5]])          # bins x seeds
    g_s = (W[:, None]*(seeds-np.array([E[b][t] for b in B[2:5]])[:, None])).sum(0)
    gap[t] = g_s.mean(); sem[t] = g_s.std(ddof=1)/np.sqrt(len(g_s))
for t in TERMS + ['mean']:
    print(f'  {t:5s} weighted gap {gap[t]:+6.0f} m (2 sem {2*sem[t]:4.0f})')
cand = {'(a) L-T coupling [C]': (gap['C'], 2*sem['C']), '(b) sediment column [sedc + dens]': (gap['sedc']+gap['dens'], 2*np.hypot(sem['sedc'], sem['dens'])),
        '(c) dyn placement [dyn]': (gap['dyn'], 2*sem['dyn'])}
order = sorted(cand, key=lambda k: cand[k][0])
print('  order (lowest-holding first; "no rank" if within 2 sem of zero or raising): ' + '; '.join(f'{k} {cand[k][0]:+.0f}' + ('' if cand[k][0] < -cand[k][1] else ' no rank') for k in order))
