"""D's planet tests (journal 66, section 5: P-a ... P-e, and the record read for P-g): the dynamic topography and
the surface flow over the last 100 Myr of each run, six-seed means per combination, against the Earth slab proxy
(SMEAN2 fast part below 300 km through the same kernels: planet_reviews.mantle_flow_earth).

  python -m planet_reviews.dynamic_topography <sweep or combination or seed folders> [--json out.json]

Bands (fixed 31 Oct, before the runs):
  P-a  share of the field's variance in degrees 1-3 >= 0.6, and P30 / P2 <= 0.1
  P-b  lows over slabs: area-weighted correlation of the field with the distance to the record's parcels younger
       than 100 Myr (below the cut) > 0
  P-c  rate at fixed points (mantle frame): median |dh/dt| <= 20 m/Myr (Rubey et al. 2017)
  P-d  rms of the air-loaded field within 0.5-2x the Earth proxy (110-460 m); its 10-Myr series over the last
       100 Myr has no step > 20%
  P-e  surface flow rms within 0.5-2x the proxy (0.7-2.7 cm/yr)
Runs with dyn_rule 'slabs' are scored on the flow's own coefficients (air-loaded, no swells: the proxy's
counterpart); the applied field (water-loaded, with swells) is reported beside them. Runs with the random field (the
control) are read from dynamic_topography_m by quadrature and are not scored for P-b / P-e.

Runs with dyn_rule 'measured' (journal 80) are checked against Earth's observed spectrum itself (OBSERVED: Valentine &
Davies 2020's coefficients, expected posterior power to l 30, converted to air loading), on dynamic_topography_air_m
(air-loaded, with the swells) by quadrature. Bands (fixed 1 Nov, before the runs; they verify the build):
  M-a  share of the power at l 1-3 within 0.5-1.5 x Earth's 0.287 (0.14-0.43)
  M-b  rms within 0.8-1.25 x Earth's 405 m (324-507 m)
  M-c  P30 / P2 <= 0.1 (Earth 0.035), as the ratio of the pooled mean powers (P2 has 5 degrees of freedom per field)
  M-d  its 10-Myr rms series has no step > 20%
  M-e  rate at fixed points, median |dh/dt| <= 20 m/Myr (Rubey et al. 2017), from snapshots <= 5 Myr apart
"""
import argparse
import json
from pathlib import Path

import numpy as np

from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_budget import seed_dirs
from planet_tectonics import mantle_flow as mf

YR_S = 3.15576e7
EARTH = {'rms_m': 222., 'share_l1_3': .84, 'p30_p2': 6e-4, 'flow_rms_cm_yr': 1.35}   # SMEAN2 slab proxy (E1)
OBSERVED = {'rms_m': 405., 'share_l1_3': .287, 'p30_p2': .035}    # Valentine & Davies 2020, air-loaded (journal 80)
M_BANDS = {'M-a share l1-3': (.14, .43), 'M-b rms m': (324., 507.), 'M-c P30/P2': (None, .1),
           'M-d max 10-Myr step': (None, .2), 'M-e rate p50 m/Myr': (None, 20.)}
BANDS = {'P-a share l1-3': (.6, None), 'P-a P30/P2': (None, .1), 'P-b corr(h, distance)': (0., None),
         'P-c rate p50 m/Myr': (None, 20.), 'P-d rms m': (110., 460.), 'P-d max 10-Myr step': (None, .2),
         'P-e flow rms cm/yr': (.7, 2.7)}
WINDOW_MYR = 100.
LMAX = 31


def _times(seed):
    return sorted(float(p.name[1:8]) for p in (Path(seed)/'snapshots').glob('t*.npz') if not p.name.endswith('-plates.json'))


def _load(seed, t):
    return np.load(Path(seed)/'snapshots'/f't{t:+07.1f}.npz')


def quadrature_alm(f, g, lmax=LMAX):
    """Coefficients of a node field by area-weighted quadrature (the cube-sphere cells as weights)."""
    return mf.point_alm(g.xyz, f*g.area/g.radius**2, lmax)


def _ok(v, band):
    lo, hi = band
    if lo is not None and not (v > lo if lo == 0. else v >= lo):     # the correlation must be strictly positive
        return False
    return hi is None or v <= hi


def seed_metrics(seed, g, coarse):
    ts = [t for t in _times(seed) if t >= -WINDOW_MYR-1e-6]
    if not ts:
        return None
    z0 = _load(seed, ts[-1])
    slabs = 'flow_topography_alm_m' in z0.files
    measured = not slabs and 'dynamic_topography_air_m' in z0.files
    tens = [t for t in ts if abs(t/10-round(t/10)) < 1e-6]
    rms, share, p30, flow, corr, rms_app, share_app, pw2, pw30 = [], [], [], [], [], [], [], [], []
    for t in tens:
        z = _load(seed, t)
        A = quadrature_alm(z['dynamic_topography_m'].astype(np.float64), g)
        P = mf.spectrum(A, LMAX)
        rms_app.append(np.sqrt(P[1:].sum())); share_app.append(P[1:4].sum()/P[1:].sum())
        if measured:
            P = mf.spectrum(quadrature_alm(z['dynamic_topography_air_m'].astype(np.float64), g), LMAX)
        if slabs:
            H = z['flow_topography_alm_m'].astype(np.complex128); lm = int(z['flow_lmax'])
            V = z['flow_velocity_alm_m_yr'].astype(np.complex128)
            P = mf.spectrum(H, lm)
            flow.append(np.sqrt(mf.gradient_power(V, lm))*100)
            young = z['slab_mass_young_kg'] > 0
            if young.any():
                h = mf.synthesis(H, lm, coarse.xyz)
                src = np.unique(fields.nearest_nodes(coarse, g.xyz[young]))
                d = np.full(coarse.size, np.inf)
                for i in range(0, src.size, 512):
                    d = np.minimum(d, np.arccos(np.clip(coarse.xyz@coarse.xyz[src[i:i+512]].T, -1, 1)).min(1))
                w = coarse.area
                hm, dm = np.average(h, weights=w), np.average(d, weights=w)
                corr.append(np.average((h-hm)*(d-dm), weights=w)/np.sqrt(np.average((h-hm)**2, weights=w)*np.average((d-dm)**2, weights=w)))
        tot = P[1:].sum()
        rms.append(np.sqrt(tot))
        share.append(P[1:4].sum()/tot if tot > 0 else np.nan); p30.append(P[30]/P[2] if P[2] > 0 else np.nan)
        pw2.append(P[2]); pw30.append(P[30])
    # P-c: consecutive 1-Myr snapshots (mantle frame: the same nodes)
    rates, rates_air = [], []
    for a, b in zip(ts[:-1], ts[1:]):
        if abs(b-a-1) > 1e-6 and not (measured and b-a <= 5+1e-6):
            continue
        za, zb = _load(seed, a), _load(seed, b)
        rates.append(np.abs(zb['dynamic_topography_m'].astype(np.float64)-za['dynamic_topography_m'])/(b-a))
        if slabs or measured:
            rates_air.append(np.abs(zb['dynamic_topography_air_m'].astype(np.float64)-za['dynamic_topography_air_m'])/(b-a))
    w = np.tile(g.area, len(rates)) if rates else None
    def wmed(r):
        if not r:
            return float('nan')
        x = np.concatenate(r); o = np.argsort(x); c = np.cumsum(w[o]); return float(x[o][np.searchsorted(c, c[-1]/2)])
    r = np.array(rms)
    steps = np.abs(np.diff(r))[r[:-1] > 0]/r[:-1][r[:-1] > 0] if len(r) > 1 else np.zeros(1)
    steps = steps if steps.size else np.zeros(1)
    m = {'slabs': slabs, 'measured': measured, 'times': tens, 'rms_series_m': [float(x) for x in rms],
         'P-a share l1-3': float(np.nanmean(share)), 'P-a P30/P2': float(np.nanmean(p30)),
         'P-b corr(h, distance)': float(np.mean(corr)) if corr else float('nan'),
         'P-c rate p50 m/Myr': wmed(rates), 'rate_air_p50_m_myr': wmed(rates_air),
         'P-d rms m': float(np.mean(rms)), 'P-d max 10-Myr step': float(steps.max()),
         'P-e flow rms cm/yr': float(np.mean(flow)) if flow else float('nan'),
         'applied_rms_m': float(np.mean(rms_app)), 'applied_share_l1_3': float(np.mean(share_app))}
    if measured:
        m.update({'M-a share l1-3': m['P-a share l1-3'], 'M-b rms m': m['P-d rms m'], 'p2_mean': float(np.mean(pw2)),
                  'p30_mean': float(np.mean(pw30)), 'M-c P30/P2': float(np.mean(pw30)/np.mean(pw2)),
                  'M-d max 10-Myr step': m['P-d max 10-Myr step'], 'M-e rate p50 m/Myr': m['rate_air_p50_m_myr']})
    if slabs:
        zl = _load(seed, tens[-1])
        m['record_kg'] = [float(x) for x in zl['slab_record_kg']]
        m['mass_by_layer_kg'] = [float(x) for x in zl['slab_mass_by_layer_kg']]
    return m


def earth_profile():
    """Excess mass of the Earth slab proxy (SMEAN2 fast part > 300 km) per model layer (depth km, kg)."""
    from planet_reviews import mantle_flow_earth as mfe
    _, loads = mfe.density_loads('smean2', True, 300.)
    return [(z, float(a[0].real*np.sqrt(4*np.pi)*(6371e3-z*1e3)**2)) for z, a in loads]


def profile_text(mass_by_layer, z_cut_km=300.):
    edges = [300, 660, 1000, 1500, 2000, 2500, 2850]
    zc = z_cut_km+25+50*np.arange(len(mass_by_layer))
    return ' '.join(f'{a}-{b}: {np.sum(np.asarray(mass_by_layer)[(zc >= a) & (zc < b)])/1e21:.2f}' for a, b in zip(edges[:-1], edges[1:]))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--json')
    a = ap.parse_args(argv)
    g, coarse = CubeSphere(217), CubeSphere(48)
    edges = [300, 660, 1000, 1500, 2000, 2500, 2850]
    try:
        ep = earth_profile()
        print(f"EARTH slab proxy (SMEAN2 fast part > 300 km): rms {EARTH['rms_m']:.0f} m, share l1-3 {EARTH['share_l1_3']:.2f}, "
              f"P30/P2 {EARTH['p30_p2']:.0e}, flow {EARTH['flow_rms_cm_yr']:.2f} cm/yr; excess mass {sum(m for _, m in ep)/1e21:.2f}e21 kg; "
              'by depth (km: 1e21 kg) ' + ' '.join(f'{lo}-{hi}: {sum(m for z, m in ep if lo <= z < hi)/1e21:.2f}' for lo, hi in zip(edges[:-1], edges[1:])))
    except FileNotFoundError:                        # data/mantle's tomography is needed for D's proxy profile only
        print('EARTH slab proxy: data/mantle tomography missing (profile not printed)')
    out = {}
    for name, seeds in seed_dirs(a.dirs):
        name = f'{Path(seeds[0]).parent.parent.name}/{name}' if seeds else name     # sweeps share combo names
        rows = {}
        for sd in seeds:
            m = seed_metrics(sd, g, coarse)
            if m is not None:
                rows[Path(sd).name] = m
        if not rows:
            continue
        slabs = any(r['slabs'] for r in rows.values())
        measured = all(r.get('measured') for r in rows.values())
        print(f"\n{name}: {len(rows)} seeds, " + ('dyn_rule slabs (scored on the flow, air-loaded)' if slabs else
              "dyn_rule measured (air-loaded field with swells; scored on M-a ... M-e against Earth's observed spectrum)"
              if measured else 'random field (applied field; P-b, P-e not applicable)'))
        for s, r in rows.items():
            print(f"  {s}: share {r['P-a share l1-3']:.2f}, P30/P2 {r['P-a P30/P2']:.1e}, corr {r['P-b corr(h, distance)']:+.2f}, "
                  f"rate {r['P-c rate p50 m/Myr']:.1f} m/Myr (air {r['rate_air_p50_m_myr']:.1f}), rms {r['P-d rms m']:.0f} m "
                  f"(max step {r['P-d max 10-Myr step']:.2f}), flow {r['P-e flow rms cm/yr']:.2f} cm/yr; applied field rms "
                  f"{r['applied_rms_m']:.0f} m, share {r['applied_share_l1_3']:.2f}")
            if r['slabs']:
                rk = r['record_kg']
                print(f"      record at {r['times'][-1]:+.0f} Myr: deposited {rk[0]:.2e}, held {rk[1]:.2e}, past 2,850 km {rk[2]:.2e} kg "
                      f"(closes to {abs(rk[0]-rk[1]-rk[2])/max(rk[0], 1):.1e}); below the cut by depth (1e21 kg) "
                      + profile_text(r['mass_by_layer_kg']))
                print('      rms series (10 Myr): ' + ' '.join(f'{x:.0f}' for x in r['rms_series_m']))
        if measured:
            mm = {k: float(np.nanmean([r[k] for r in rows.values()])) for k in M_BANDS}
            mm['M-c P30/P2'] = float(np.sum([r['p30_mean'] for r in rows.values()])/np.sum([r['p2_mean'] for r in rows.values()]))
            vm = {k: 'pass' if _ok(mm[k], band) else 'FAIL' for k, band in M_BANDS.items()}
            print(f"  Earth observed (air-loaded): rms {OBSERVED['rms_m']:.0f} m, share l1-3 {OBSERVED['share_l1_3']:.3f}, "
                  f"P30/P2 {OBSERVED['p30_p2']:.3f}")
            print('  six-seed mean: ' + '; '.join(f"{k} {mm[k]:.3g} {vm[k]}" for k in M_BANDS)
                  + f"; applied field rms {np.mean([r['applied_rms_m'] for r in rows.values()]):.0f} m")
            print(f"  M-a ... M-e on the mean: {'all pass' if all(v == 'pass' for v in vm.values()) else 'not all pass'}")
            out[name] = {'seeds': rows, 'mean': mm, 'verdict': vm}
            continue
        mean = {k: float(np.nanmean([r[k] for r in rows.values()])) for k in BANDS}
        verdict = {}
        for k, band in BANDS.items():
            if not slabs and k in ('P-b corr(h, distance)', 'P-e flow rms cm/yr'):
                verdict[k] = 'n/a'; continue
            verdict[k] = 'pass' if _ok(mean[k], band) else 'FAIL'
        print('  six-seed mean: ' + '; '.join(f"{k} {mean[k]:.3g} {verdict[k]}" for k in BANDS))
        if slabs:
            allp = all(v == 'pass' for v in verdict.values())
            print(f"  P-a ... P-e on the mean: {'all pass' if allp else 'not all pass'} (adoption also needs the gates on every seed)")
        out[name] = {'seeds': rows, 'mean': mean, 'verdict': verdict}
    if a.json:
        Path(a.json).write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
