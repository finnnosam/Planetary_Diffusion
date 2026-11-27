"""Step-back probe (journal 98): split the land's height excess into the sea, the crust-thickness distribution and
the column (elevation at a given crust + sediment thickness). Six seeds of a sweep arm, at -250 Myr and 0 Myr.
Earth: CRUST1.0 x ETOPO thickness bins (earth-crust-elevation-n217.json), ETOPO land hypsometry, Seton ages.
T bins carry CRUST1.0's thickness error (regression dilution, decision 1 Nov): read the spreads beside them."""
import sys, json, numpy as np
arm = sys.argv[1] if len(sys.argv) > 1 else 'outputs/planet-calibration/sweep-20261108-124109/rift_line-force'
REF = 'outputs/planet-reference/'
B = json.load(open(REF + 'earth-crust-elevation-n217.json'))['bins']
edges = [tuple(b['T_km']) for b in B]; efrac = np.array([b['area_frac'] for b in B])
ez = [b['earth_z_p25_p50_p75_m'][1] for b in B]
E = np.load(REF + 'earth-n217.npz')['elevation_m']; age_e = np.load(REF + 'earth-seafloor-age-n217.npz')['age_myr']
el = E > 0
print('Earth land %.3f, p25/p50/p90 %s, p90-p25 %.0f' % (el.mean(), np.round(np.percentile(E[el], [25, 50, 90])),
      np.subtract(*np.percentile(E[el], [90, 25]))))
def wpct(x, w, q):
    o = np.argsort(x); c = np.cumsum(w[o]); c /= c[-1]; return np.interp(q, c, x[o])
AGES = [(2, 10), (10, 30), (30, 60), (60, 100), (100, 200)]
print('Earth floor by age (p50, m): ' + ' '.join(f'{a}-{b}: {np.nanmedian(E[(age_e >= a) & (age_e < b)]):.0f}' for a, b in AGES))
for tag in ('t-0250.0', 't+0000.0'):
    Z, T, OA, LA, SL = [], [], [], [], []
    for s in range(1001, 1007):
        d = np.load(f'{arm}/seed{s}/snapshots/{tag}.npz'); sl = float(d['sea_level_m'])
        Z.append(d['elevation_m'] + sl); T.append((d['crust_thickness_m'] + d['sediment_m']) / 1e3)
        OA.append(d['ocean_age_myr']); LA.append(d['lithosphere_age_myr']); SL.append(sl)
    Z, T, OA, LA = map(np.concatenate, (Z, T, OA, LA)); sea = float(np.mean(SL))   # Z on the datum
    w = np.zeros_like(Z)
    print(f'\n== {tag}: sea on the datum {sea:.0f} m')
    print('  T bin | area planet / Earth | z on datum p50 | z above sea p50 | Earth p50 | old lith (> 1 Gyr) z above sea')
    for (a, b), f, e in zip(edges, efrac, ez):
        m = (T >= a) & (T < b)
        if not m.any(): continue
        w[m] = f / m.mean(); mo = m & (LA > 1000)
        print(f'  {a}-{b} | {m.mean():.4f} / {f:.4f} | {np.median(Z[m]):.0f} | {np.median(Z[m]) - sea:.0f} | {e:.0f} | '
              + (f'{np.median(Z[mo]) - sea:.0f}' if mo.sum() > 50 else '-'))
    oc = T < 15
    print('  floor by age, on datum (p50): ' + ' '.join(f'{a}-{b}: {np.nanmedian(Z[oc & (OA >= a) & (OA < b)]):.0f}' for a, b in AGES))
    for name, s0, wt in (('as run', sea, None), ('sea at the datum', 0., None), ('as run, Earth T weights', sea, w),
                         ('sea at the datum, Earth T weights', 0., w)):
        l = Z > s0; ww = np.ones(l.sum()) if wt is None else wt[l]
        q = [wpct(Z[l] - s0, ww, x) for x in (.25, .5, .9)]
        lf = l.mean() if wt is None else (wt * l).sum() / wt.sum()
        print(f'  {name}: land p25/p50/p90 {np.round(q)}, p90-p25 {q[2] - q[0]:.0f}' + (f', land {lf:.3f}' if wt is None else ''))
