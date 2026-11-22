"""Borrowed control (the hundredth pass, journal 134): scripts.sweep_plates --borrow NAME=FOLDER re-runs only the
control's check seed, compares it with the earlier arm's run field by field, and on a match places the other seeds
from that arm (hard links for the data files); on a mismatch it runs them itself. The runs and Earth's scores are
stubbed: these tests check the bookkeeping, not the model."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from scripts import sweep_plates as spl
from planet_reviews.hypsometry import QS


def fake_run(seed_dir, seed, offset=0.):
    """A seed folder as a run writes it: snapshots (npz + plates json), stats.json, a log and a picture."""
    sd = Path(seed_dir)
    (sd/'snapshots').mkdir(parents=True, exist_ok=True)
    z = np.arange(12, dtype=np.float32)+seed+offset
    np.savez(sd/'snapshots'/'t+0000.0.npz', elevation_m=z, continental=z > 5)
    np.savez(sd/'snapshots'/'t-0005.0.npz', elevation_m=z-1, continental=z > 6)
    (sd/'snapshots'/'t+0000.0-plates.json').write_text(json.dumps({'plates': [1, 2]}))
    (sd/'stats.json').write_text(json.dumps({'seed': seed, 'years': 250., 'seconds': 1.,
                                             'params_effective': {'flow_power': 'measured'}}))
    (sd/'run.log').write_text('log\n')
    (sd/'present.png').write_bytes(b'png')


def test_plan_runs_only_the_check_seed_of_the_borrowed_combination_and_first():
    combos = [{'k': 'a'}, {'k': 'b'}]
    names = ['k-a', 'k-b']
    plan = spl.plan_jobs(combos, names, [1001, 1002, 1003], borrow='k-a')
    assert plan[0] == (combos[0], 1001)
    assert [s for c, s in plan if c is combos[0]] == [1001]
    assert [s for c, s in plan if c is combos[1]] == [1001, 1002, 1003]
    assert spl.plan_jobs(combos, names, [1001, 1002]) == [(c, s) for c in combos for s in (1001, 1002)]
    two = spl.plan_jobs(combos, names, [1001, 1002, 1003], borrow='k-b', check_seeds=[1002, 1003])
    assert [s for c, s in two if c is combos[1]] == [1002, 1003] and two[0] == (combos[1], 1002)


def test_borrow_seed_links_the_data_and_copies_the_rest(tmp_path):
    src, dst = tmp_path/'ref'/'seed1002', tmp_path/'new'/'seed1002'
    fake_run(src, 1002)
    nl, nc = spl.borrow_seed(src, dst)
    files = sorted(p.relative_to(src) for p in src.rglob('*') if p.is_file())
    assert files == sorted(p.relative_to(dst) for p in dst.rglob('*') if p.is_file())
    assert nl == 4 and nc == 2                       # two npz, the plates json, stats.json | run.log, present.png
    for f in files:
        assert (src/f).read_bytes() == (dst/f).read_bytes()
        assert os.path.samefile(src/f, dst/f) == (f.suffix in spl.LINKED_SUFFIXES)


def test_same_run_reads_compare_runs_verdict(tmp_path):
    fake_run(tmp_path/'a', 1001)
    fake_run(tmp_path/'b', 1001)
    fake_run(tmp_path/'c', 1001, offset=0.5)
    ok, txt = spl.same_run(tmp_path/'a', tmp_path/'b')
    assert ok and '0 differ' in txt
    ok, txt = spl.same_run(tmp_path/'a', tmp_path/'c')
    assert not ok and 'elevation_m' in txt


def test_a_source_missing_a_seed_or_run_differently_is_refused(tmp_path):
    arm = tmp_path/'sweep-old'/'arm'
    for s in (1001, 1002):
        fake_run(arm/f'seed{s}', s)
    (tmp_path/'sweep-old'/'sweep.json').write_text(json.dumps({'years': 250., 'save_every': 5}))
    spl.check_borrow_source(arm, [1001, 1002], 250., 5)
    with pytest.raises(SystemExit):
        spl.check_borrow_source(arm, [1001, 1002, 1003], 250., 5)
    with pytest.raises(SystemExit):
        spl.check_borrow_source(arm, [1001, 1002], 250., 1)


class _Tracker:
    def __init__(self):
        self.lines = []

    def log(self, s):
        self.lines.append(s)

    def update(self, *a, **k):
        pass


def _stub(monkeypatch, tmp_path, drift_seed=None):
    """sweep_plates.main without the model or Earth: each job writes fake_run (offset on drift_seed's control run)."""
    import planet_grid.sphere_graph as sg
    import planet_reviews.plate_texture as pt
    import planet_earth.reference as ref
    from scripts import progress as pg
    monkeypatch.setattr(spl, 'OUTPUTS', tmp_path/'outputs')
    monkeypatch.setattr(sg, 'CubeSphere', lambda n: None)
    monkeypatch.setattr(pt, 'earth_texture', lambda g: {'band_power_m2': {'b': 1.}, 'relief_368km_m': {'25': 1., '50': 1., '75': 1.}})
    monkeypatch.setattr(ref, 'load', lambda n: {'no_ice_sheets': {'land_fraction': .28}})
    monkeypatch.setattr(pg, 'record', lambda *a, **k: None)
    ran = []

    def run_jobs(jobs, workers, years, tracker, label):
        for i, (c, s, out, yrs, se) in enumerate(jobs):
            ran.append((c['flow_power'], s))
            off = 0.5 if (s == drift_seed and c['flow_power'] == 'flow') else 0.
            fake_run(Path(out)/f'seed{s}', s, off)
            yield i, {'params': c, 'seed': s, 'dir': str(Path(out)/f'seed{s}'), 'seconds': 60.}
    monkeypatch.setattr(spl, 'run_jobs', run_jobs)
    shelf = {'whole_p70_m': -40., 'continental_submerged': .3, 'land_p25_m': 200., 'shallow_sea_share': .05,
             'continental_fraction': .4}
    monkeypatch.setattr(spl, 'score', lambda d, et, eh, g: {
        'score': .3, 'texture_score': .2, 'hypsometry_score': .1, 'arrangement_score': .1, 'budget_score': .01,
        'flux_score': 0., 'export_km3_yr': 6., 'continental_volume_change': .0,
        'arrangement_ratios': {k: 1. for k in spl.ARRANGEMENT_KEYS},
        'hypsometry': {'land_fraction': .28, **{q: 100. for q in QS}},
        'texture': {'band_power_ratio': {'b': 1.}, 'relief_ratio': {'25': 1., '50': 1., '75': 1.}},
        'shelf': shelf, 'shelf_pass': True})
    ref_arm = tmp_path/'outputs'/'planet-calibration'/'sweep-old'/'flow_power-flow'
    for s in (1001, 1002, 1003):
        fake_run(ref_arm/f'seed{s}', s)
    (ref_arm.parent/'sweep.json').write_text(json.dumps({'years': 250., 'save_every': 5}))
    return ref_arm, ran


@pytest.mark.parametrize('drift', [False, True])
def test_sweep_borrows_on_a_match_and_runs_the_control_on_a_mismatch(monkeypatch, tmp_path, drift):
    ref_arm, ran = _stub(monkeypatch, tmp_path, drift_seed=1001 if drift else None)
    tr = _Tracker()
    out = spl.main(['--grid', 'flow_power=flow,measured', '--seeds', '1001,1002,1003', '--rule', 'shelf',
                    '--borrow', f'flow_power-flow={ref_arm}'], tracker=tr)
    ctl = out/'flow_power-flow'
    info = json.loads((ctl/'borrowed.json').read_text())
    meta = json.loads((out/'sweep.json').read_text())['borrow']
    res = json.loads((out/'results.json').read_text())
    assert ran[0] == ('flow', 1001)                  # the check seed is scheduled first
    assert sorted(s for c, s in ran if c == 'measured') == [1001, 1002, 1003]
    assert sorted(r['seed'] for r in res['runs'] if r['params']['flow_power'] == 'flow') == [1001, 1002, 1003]
    assert {r['seed'] for r in res['runs']} == {1001, 1002, 1003} and len(res['runs']) == 6
    assert len(res['table']) == 2
    if not drift:
        assert info['identical'] and info['borrowed_seeds'] == [1002, 1003] and meta['identical']
        assert [s for c, s in ran if c == 'flow'] == [1001]
        assert os.path.samefile(ctl/'seed1002'/'stats.json', ref_arm/'seed1002'/'stats.json')
        assert 'IDENTICAL' in (ctl/'borrowed.txt').read_text()
    else:
        assert not info['identical'] and info['borrowed_seeds'] == [] and not meta['identical']
        assert [s for c, s in ran if c == 'flow'] == [1001, 1002, 1003]
        assert not os.path.samefile(ctl/'seed1002'/'stats.json', ref_arm/'seed1002'/'stats.json')
        assert 'DIFFERENT' in (ctl/'borrowed.txt').read_text()
        assert any('DIFFER' in l for l in tr.lines)
