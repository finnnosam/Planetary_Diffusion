"""Pass 84 (journal 109): the thick tail by arm. Label-free continent (crust + sediment >= 20 km, as crust_review)
by thickness band, at 0 Myr and as the mean over the saved snapshots -40 ... 0 Myr; Earth on CRUST1.0's 1 degree
cells (cos-lat weights). Beside it: the steady Moho temperature of the planet's continental columns (crust_flow's
geotherm on the crust and C's L), the elevation of crust > 45 km, and the crust flow's book (substeps per Myr, share
of edges at the numerical ceiling) where the run used the measured rheology.
Usage (repo root): python thick_tail.py <arm> [<arm> ...] [--seeds 1001,...,1006]"""
import json
import sys
from pathlib import Path

import numpy as np
from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_reviews.crust_review import column_thickness_m, label_free_continent, run_params, wpct
from planet_tectonics.crust_flow import moho_state
from planet_tectonics.lithosphere import frame, thickness_km

g = CubeSphere(217)
BANDS = [(20, 30), (30, 35), (35, 40), (40, 45), (45, 50), (50, 55), (55, 65), (65, 1e9)]
TIMES = ['t-0040.0', 't-0035.0', 't-0030.0', 't-0025.0', 't-0020.0', 't-0015.0', 't-0010.0', 't-0005.0', 't+0000.0']


def dist(T_km, w):
    w = w/w.sum()
    return [float(w[(T_km >= lo) & (T_km < hi)].sum()) for lo, hi in BANDS]


def earth():
    T = np.loadtxt(DATA/'crust1'/'crsthk.xyz')[:, 2].reshape(180, 360).ravel()*1e3
    a = np.repeat(np.cos(np.radians(89.5-np.arange(180)))[:, None], 360, 1).ravel()
    c = label_free_continent(T)
    return dist(T[c]/1e3, a[c])


def snap(seed_dir, t):
    f = Path(seed_dir)/'snapshots'/f'{t}.npz'
    if not f.exists():
        return None
    z = np.load(f)
    p = run_params(seed_dir)
    Tcol = column_thickness_m(z, p)
    c = label_free_continent(Tcol)
    out = {'dist': dist(Tcol[c]/1e3, g.area[c]), 'cont_area': float(g.area[c].sum()/g.area.sum())}
    e = z['elevation_m'].astype(float)
    thick = c & (Tcol > 45e3)
    out['z_thick_p50'] = wpct(e[thick], g.area[thick], 50) if thick.any() else np.nan
    flag = z['continental'].astype(bool)
    Tc = z['crust_thickness_m'].astype(float)
    if 'lithosphere_age_myr' in z.files and p.get('lithosphere') == 'thermal_age':
        age = z['lithosphere_age_myr'].astype(float)
        L = np.where(np.isfinite(age), thickness_km(np.nan_to_num(age, nan=0.), p), frame(p)[0])*1e3
    else:
        L = np.full(g.size, frame(p)[0]*1e3)
    fp = {'flow_heat_upper_w_m3': 1.0e-6, 'flow_heat_lower_w_m3': 0.4e-6, 'flow_k_crust': 2.5, 'flow_k_mantle': 3.0,
          'flow_t_surface_c': 10., 'flow_t_lab_c': 1175., **p}
    _, _, tm = moho_state(Tc[flag], L[flag], fp)
    w = g.area[flag]/g.area[flag].sum()
    out['moho_gt700'] = float(w[tm > 700].sum()); out['moho_gt800'] = float(w[tm > 800].sum())
    out['L_thick_p50'] = wpct(L[flag & (Tc > 45e3)]/1e3, g.area[flag & (Tc > 45e3)], 50) if (flag & (Tc > 45e3)).any() else np.nan
    return out


def book(seed_dir):
    st = json.loads((Path(seed_dir)/'stats.json').read_text())
    ser = st['series']
    e0, e1 = ser[0]['erosion_totals_km3'], ser[-1]['erosion_totals_km3']
    span = ser[-1]['time_myr']-ser[0]['time_myr']
    sub = (e1.get('crust_flow_substeps', 0.)-e0.get('crust_flow_substeps', 0.))*1e9
    ed = (e1.get('crust_flow_edges', 0.)-e0.get('crust_flow_edges', 0.))
    cap = (e1.get('crust_flow_capped_edges', 0.)-e0.get('crust_flow_capped_edges', 0.))
    return sub/span if span else np.nan, cap/ed if ed else np.nan


def main(argv):
    seeds = list(range(1001, 1007))
    if '--seeds' in argv:
        i = argv.index('--seeds'); seeds = [int(s) for s in argv[i+1].split(',')]; argv = argv[:i]+argv[i+2:]
    hdr = ' / '.join(f'{lo:.0f}-{hi:.0f}' if hi < 1e8 else f'>{lo:.0f}' for lo, hi in BANDS)
    E = earth()
    print(f'label-free continent by crust + sediment (km): {hdr}   | > 45 / > 55')
    print(f'  Earth (CRUST1.0, 1 deg)        : ' + ' / '.join(f'{v:.3f}' for v in E) + f'   | {sum(E[4:]):.3f} / {sum(E[6:]):.3f}')
    for arm in argv:
        name = Path(arm).name
        for lab, times in (('0 Myr', TIMES[-1:]), ('window', TIMES)):
            rs = [r for s in seeds for t in times if (r := snap(Path(arm)/f'seed{s}', t)) is not None]
            if not rs:
                continue
            D = np.mean([r['dist'] for r in rs], 0)
            sd = np.std([sum(r['dist'][4:]) for r in rs if True]) if lab == '0 Myr' else np.nan
            print(f'  {name[:22]:22s} {lab:7s}: ' + ' / '.join(f'{v:.3f}' for v in D) +
                  f'   | {D[4:].sum():.3f} / {D[6:].sum():.3f}' + (f' (sd > 45 {sd:.3f})' if lab == '0 Myr' else ''))
            print(f'      area {np.mean([r["cont_area"] for r in rs]):.3f}; elevation of crust > 45 km p50 '
                  f'{np.nanmean([r["z_thick_p50"] for r in rs]):.0f} m; L under crust > 45 km p50 '
                  f'{np.nanmean([r["L_thick_p50"] for r in rs]):.0f} km; steady Moho > 700 / > 800 C '
                  f'{np.mean([r["moho_gt700"] for r in rs]):.3f} / {np.mean([r["moho_gt800"] for r in rs]):.3f} of the continent')
        bk = [book(Path(arm)/f'seed{s}') for s in seeds if (Path(arm)/f'seed{s}'/'stats.json').exists()]
        if bk and np.isfinite(bk[0][1]):
            print(f'      crust flow: substeps per Myr {np.mean([b[0] for b in bk]):.0f}, edges at the ceiling '
                  f'{np.mean([b[1] for b in bk]):.4f}')


if __name__ == '__main__':
    main(sys.argv[1:])
