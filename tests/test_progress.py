"""Sweep progress and timing record (sixty-eighth pass, journal 90): output only."""
import io
import json

import numpy as np

from scripts import progress as pg


def test_fmt_and_schedule():
    assert pg.fmt(0) == '0:00:00' and pg.fmt(3725.4) == '1:02:05'
    # 24 runs of 20 min on 8 workers: three waves
    assert pg.schedule_left([], 24, 1200., 8) == 3600.
    # 12 runs: two waves, the second of four
    assert pg.schedule_left([], 12, 1200., 8) == 2400.
    # running runs first, then the queue onto the worker free first
    assert pg.schedule_left([100., 600.], 1, 1000., 2) == 1100.


def test_record_and_predict(tmp_path):
    p = tmp_path/'t.jsonl'
    assert pg.predict('pytest', path=p) == pg.DEFAULTS[('pytest', None)]
    for s in (100., 110., 130.):
        pg.record('pytest', 'x', s, path=p, files=19)
    assert pg.predict('pytest', path=p) == 110.
    pg.record('review', 'coast_budget', 24., path=p, runs=12)
    assert pg.predict('review', 'coast_budget', per='runs', path=p) == 2.
    assert pg.predict_run(8, 250., path=p) == pg.DEFAULTS[('run', None)]
    assert pg.predict_run(8, 125., path=p) == pg.DEFAULTS[('run', None)]/2
    pg.record('run', 'a', 900., path=p, workers=8, years=250.)
    pg.record('run', 'a', 1100., path=p, workers=8, years=250.)
    assert pg.predict_run(8, 250., path=p) == 1000.
    assert pg.predict_run(4, 125., path=p) == 500.        # other setups: scaled by length
    with open(p, 'a') as f:
        f.write('not json\n')
    assert len(pg.read_records(p)) == 6


def test_tracker_plain_output():
    s = io.StringIO()
    tr = pg.Tracker([('tests', 10.), ('sweep', 100.)], stream=s, live=False, width=200)
    tr.start_step()
    tr.update(.5, '5 of 10 files')
    tr.log('tests passed')
    tr.end_step()
    tr.start_step()
    tr.update(.25, '1 of 4 runs done')
    tr.end_step()
    tr.finish('done')
    out = s.getvalue()
    assert 'tests passed\n' in out and out.rstrip().endswith('done')
    assert '1/2 tests' in out and '2/2 sweep: 1 of 4 runs done' in out and '\r' not in out
    assert len(tr.summary()) == 4 and tr.summary()[-1].startswith('total')


def test_tracker_live_redraws_in_place():
    s = io.StringIO()
    tr = pg.Tracker([('sweep', 60.)], stream=s, live=True, width=120)
    tr.start_step()
    tr.update(.5, 'half')
    tr.log('run 1 of 2 done')
    out = s.getvalue()
    assert out.startswith('\r[') and 'run 1 of 2 done\n' in out
    assert all(len(l) <= 120 for l in out.replace('\r', '\n').split('\n'))


def test_review_commands():
    from pathlib import Path
    from scripts import run_sweep as rs
    sw = Path('outputs')/'planet-calibration'/'sweep-x'
    cmds = {r[0]: rs.review_cmd(r, sw, 'arm-a', None) for r in rs.REVIEWS}
    assert cmds['identity'] is None                       # no reference given
    assert str(sw/'arm-a') in cmds['coast_budget'] and cmds['coast_budget'][-2] == '--json'
    assert cmds['plate_pieces_review'][-1] == 'arm-a'
    assert rs.review_cmd(rs.REVIEWS[4], sw, None, None) is None   # coast_budget needs the control
    assert all(Path(f'tests/{t}.py').exists() for t in rs.PLANET_TESTS)


def test_progress_hook_is_output_only(tmp_path):
    from scripts import simulate_plates as sp
    seen = []
    argv = ['--n', '48', '--allow-coarse', '--seeds', '3', '--years', '3', '--param', 'threads=1']
    sp.main(argv+['--out', str(tmp_path/'a')])
    sp.PROGRESS = seen.append
    try:
        sp.main(argv+['--out', str(tmp_path/'b')])
    finally:
        sp.PROGRESS = None
    assert seen == [-3., -2., -1., 0.]
    a = json.loads((tmp_path/'a'/'seed3'/'stats.json').read_text())
    b = json.loads((tmp_path/'b'/'seed3'/'stats.json').read_text())
    a.pop('seconds'); b.pop('seconds')
    assert a == b
    za = np.load(tmp_path/'a'/'seed3'/'snapshots'/'t+0000.0.npz')
    zb = np.load(tmp_path/'b'/'seed3'/'snapshots'/'t+0000.0.npz')
    assert all(np.array_equal(za[k], zb[k], equal_nan=True) for k in za.files)
