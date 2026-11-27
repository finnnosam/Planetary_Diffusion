"""Pass 105 (journal 142): the orogeny-age hand-off replay, read. One saved plate history (the review's Chain run, seed
1001) through the landscape and the meso twice on the same forcing: the legacy hand-off (active belts eroded as old rock,
the review's defect) and the plate layer's own law (`--age-rule plate`, the fix). The legacy replay is also compared
with the Chain run's own landscape (same code otherwise): identical means the landscape repeats and the difference
between the two replays is the age rule alone.

Measures (fixed in journal 142 before the runs): erodibility on land above 3 km (mean; share >= 0.5); land fraction,
land p50 / p90 / p99, land above 3 / 4 km and the largest connected area above 3 / 4 km (eight-neighbour graph of the
n = 868 grid); denudation by elevation (last 10 Myr); rivers: river-node length (nodes x spacing), median and p99
drainage area of river nodes, land draining to the sea; 368 km windows: class shares (flat, plain-upland, hills,
mountains, plains-lowlands), the input's curvature per class; meso detail given the class: incision and 5-11 km band
ratios against Earth (pass >= 0.7, over-textured > 2).

  python replay_compare.py <chain run folder> <replay folder> > <replay folder>\\compare.txt
"""
import json
import sys
from pathlib import Path

import numpy as np

from planet_grid.sphere_graph import CubeSphere

CLASSES = ('flat', 'plain-upland', 'hills', 'mountains', 'plains-lowlands')


def wq(v, w, qs):
    o = np.argsort(v); c = np.cumsum(w[o]); c /= c[-1]
    return [float(v[o][min(np.searchsorted(c, q/100), len(v)-1)]) for q in qs]


def largest_component(g, m):
    """Largest connected area (km^2) of mask m on the graph's neighbours (label propagation on the masked nodes)."""
    idx = np.flatnonzero(m)
    if idx.size == 0:
        return 0.
    nb = g.nbr[idx]
    local = np.full(g.size, -1, np.int64); local[idx] = np.arange(idx.size)
    nl = np.where(nb >= 0, local[np.clip(nb, 0, None)], -1)
    lab = np.arange(idx.size)
    while True:
        cand = np.where(nl >= 0, lab[np.clip(nl, 0, None)], idx.size)
        new = np.minimum(lab, cand.min(1))
        new = new[new]                                   # pointer jumping
        if np.array_equal(new, lab):
            break
        lab = new
    area = np.bincount(lab, weights=g.area[idx])
    return float(area.max()/1e6)


def landscape_read(d, g):
    z = np.load(d/'present.npz')
    e = z['elevation_m'].astype(np.float64); oc = z['ocean']; land = ~oc
    a = g.area
    s = json.loads((d/'stats.json').read_text())
    er = z['erodibility'].astype(np.float64)
    hi3, hi4 = land & (e > 3000), land & (e > 4000)
    riv = z['river'] & land
    A = z['drainage_area_m2'].astype(np.float64)
    spacing_km = np.pi*g.radius/2/g.n/1e3
    r = {'age_rule': s.get('age_rule', 'legacy (before journal 142)'),
         'land': float(a[land].sum()/a.sum()),
         'land_p50_p90_p99': wq(e[land], a[land], (50, 90, 99)),
         'land_gt3km_of_land': float(a[hi3].sum()/a[land].sum()), 'land_gt4km_of_land': float(a[hi4].sum()/a[land].sum()),
         'largest_gt3km_km2': largest_component(g, hi3), 'largest_gt4km_km2': largest_component(g, hi4),
         'erod_land_gt3km_mean': float(er[hi3].mean()) if hi3.any() else float('nan'),
         'erod_land_gt3km_share_ge05': float((er[hi3] >= .5).mean()) if hi3.any() else float('nan'),
         'erod_land_mean': float(er[land].mean()),
         'denudation': s['denudation_last10myr_m_per_myr'], 'denudation_score': s['denudation_score'],
         'river_length_km': float(riv.sum()*spacing_km),
         'river_area_p50_p99_km2': [float(np.percentile(A[riv], q)/1e6) for q in (50, 99)] if riv.any() else [float('nan')]*2,
         'land_draining_to_sea': s.get('land_draining_to_sea'),
         'continental_crust_submerged': s.get('continental_crust_submerged')}
    return r, e


def meso_read(d):
    if not (d/'results.json').exists():
        return None
    r = json.loads((d/'results.json').read_text())
    inp = r['inputs']['planet']
    return {'shares': inp['class_share_of_land_windows'], 'curvature': inp['curvature_median_m'],
            'verdict': {k: {kk: v.get(kk) for kk in ('incision_ratio', 'band_5_11_ratio', 'passes', 'over_textured')}
                        for k, v in r['verdict'].items()},
            'earth_shares': r['inputs']['earth_meso_val']['class_share_of_land_windows']}


def main(chain, replay):
    chain, replay = Path(chain), Path(replay)
    g = CubeSphere(868)
    arms = {'chain (review, legacy)': (chain/'landscape', chain/'meso'),
            'replay legacy': (replay/'legacy'/'landscape', replay/'legacy'/'meso'),
            'replay plate (fix)': (replay/'plate'/'landscape', replay/'plate'/'meso')}
    L, M, Z = {}, {}, {}
    for k, (ld, md) in arms.items():
        if (ld/'present.npz').exists():
            L[k], Z[k] = landscape_read(ld, g)
            M[k] = meso_read(md)
    print(f'chain {chain}; replay {replay}')
    if 'chain (review, legacy)' in Z and 'replay legacy' in Z:
        dz = np.abs(Z['chain (review, legacy)']-Z['replay legacy'])
        print(f'repeat check (legacy replay against the Chain run): max |dz| {dz.max():.3g} m, '
              f'{(dz > 0).mean():.3%} of nodes differ' + ('  IDENTICAL' if dz.max() == 0 else ''))
    keys = list(L)
    print('\nlandscape (n = 868) | ' + ' | '.join(keys))
    rows = [('erodibility on land > 3 km, mean', lambda r: f"{r['erod_land_gt3km_mean']:.3f}"),
            ('erodibility on land > 3 km, share >= 0.5', lambda r: f"{r['erod_land_gt3km_share_ge05']:.3f}"),
            ('erodibility on land, mean', lambda r: f"{r['erod_land_mean']:.3f}"),
            ('land fraction', lambda r: f"{r['land']:.4f}"),
            ('land p50 / p90 / p99 (m)', lambda r: ' / '.join(f'{v:.0f}' for v in r['land_p50_p90_p99'])),
            ('land > 3 km / > 4 km, of land', lambda r: f"{r['land_gt3km_of_land']:.4f} / {r['land_gt4km_of_land']:.4f}"),
            ('largest area > 3 km / > 4 km (Mkm2)', lambda r: f"{r['largest_gt3km_km2']/1e6:.3f} / {r['largest_gt4km_km2']/1e6:.3f}"),
            ('denudation mean land / 200-500 / 500-1000 / 1000-2000 (m/Myr)',
             lambda r: ' / '.join(f"{r['denudation'][k]:.0f}" for k in ('mean_land', '200-500', '500-1000', '1000-2000'))),
            ('denudation score', lambda r: f"{r['denudation_score']:.2f}"),
            ('river-node length (km)', lambda r: f"{r['river_length_km']:.4g}"),
            ('river drainage area p50 / p99 (km2)', lambda r: ' / '.join(f'{v:.4g}' for v in r['river_area_p50_p99_km2'])),
            ('land draining to the sea', lambda r: f"{r['land_draining_to_sea']:.3f}"),
            ('continental crust submerged', lambda r: f"{r['continental_crust_submerged']:.3f}")]
    for name, fn in rows:
        print(f'{name} | ' + ' | '.join(fn(L[k]) for k in keys))
    mk = [k for k in keys if M.get(k)]
    if mk:
        e = M[mk[0]]['earth_shares']
        print('\n368 km land windows, class share | ' + ' | '.join(mk) + ' | Earth (meso validation)')
        for c in CLASSES:
            print(f'{c} | ' + ' | '.join(f"{M[k]['shares'][c]:.4f}" for k in mk) + f' | {e[c]:.4f}')
        print('\ninput curvature median (m) | ' + ' | '.join(mk))
        for c in CLASSES:
            print(f'{c} | ' + ' | '.join(f"{M[k]['curvature'].get(c, float('nan')):.1f}" for k in mk))
        print('\nmeso detail given the class: incision ratio / 5-11 km band ratio (pass >= 0.7; over > 2) | ' + ' | '.join(mk))
        for c in CLASSES:
            print(f'{c} | ' + ' | '.join(
                f"{M[k]['verdict'][c]['incision_ratio']:.2f} / {M[k]['verdict'][c]['band_5_11_ratio']:.2f}"
                f"{' pass' if M[k]['verdict'][c]['passes'] else ' fail'}{' over' if M[k]['verdict'][c]['over_textured'] else ''}"
                for k in mk))
    (replay/'compare.json').write_text(json.dumps({'landscape': L, 'meso': M}, indent=1))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
