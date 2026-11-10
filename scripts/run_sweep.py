"""A whole sweep under one progress line: the tests, the runs, the reviews (sixty-eighth pass, journal 90; the user
asked on 7 Nov for elapsed against predicted time, overall and step percentages and plain 'run k done' lines).

  python -m scripts.run_sweep --name craton-lithosphere --grid lithosphere=none,thermal_age \\
      --control lithosphere-none --identity outputs\\planet-calibration\\sweep-20261107-094624\\<arm>

Steps: 1 tests (pytest, the planet-layer subset PLANET_TESTS unless --tests is given; its printout in
outputs/planet-calibration/<name>-pytest.txt and <sweep>/reviews/pytest.txt; a failure stops here and prints it);
2 sweep (scripts.sweep_plates with --grid / --seeds / --rule / --save-every / --workers); 3 reviews (REVIEWS, each
into <sweep>/reviews/<review>.txt, and .json where it writes one). The console shows one line, redrawn in place:

  [#########-----------]  45% | 0:27:10 of ~1:00:20 | 2/3 sweep: 9 of 24 runs done, 8 running [41%]

and a line per finished run or review. The prediction comes from outputs/planet-calibration/timings.jsonl (every
step's, run's and review's time, written as they finish; scripts.progress), so it follows the machine. At the end
<sweep>/reviews/timings.txt lists each step's predicted and actual time. Exit status 1 on any failure.
"""
import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

from planet_grid.paths import OUTPUTS
from scripts import progress as pg

PLANET_TESTS = ('test_shelf', 'test_plates', 'test_sphere_graph', 'test_graph_landscape', 'test_coupling',
                'test_boundary_lengths', 'test_erosion_law_check', 'test_crust_reviews', 'test_trench_pairs',
                'test_mantle_flow', 'test_flow_reference', 'test_save_cadence', 'test_lithosphere', 'test_dyn_measured',
                'test_margin_wedge', 'test_collision_spread', 'test_breakup_match', 'test_start_breakup',
                'test_init_site', 'test_progress', 'test_rift_speed', 'test_rift_line', 'test_pelagic', 'test_deep_crust',
                'test_collision_lithosphere', 'test_lith_frame')
# (name, module or script, arguments; {SW} sweep folder, {CTL} control folder, {CTLNAME} its name, {REF} identity
# reference, {OUT} the review's json path): as the sweep .cmd files ran them since pass 67
REVIEWS = (
    ('identity', 'tools/compare_runs.py', ['{CTL}', '{REF}']),
    ('max_elevation', 'tools/max_elevation.py', ['{SW}']),
    ('rift_production', 'planet_reviews.rift_production', ['{SW}', '--json', '{OUT}']),
    ('passive_survival', 'planet_reviews.passive_survival', ['{SW}', '--json', '{OUT}']),
    ('coast_budget', 'planet_reviews.coast_budget', ['{SW}', '{CTL}', '--json', '{OUT}']),
    ('coast_profile', 'planet_reviews.coast_profile', ['{SW}', '{CTL}', '--json', '{OUT}']),
    ('plate_pieces_review', 'scripts.plate_pieces_review', ['--latest', '--control', '{CTLNAME}']),
    ('crust_review', 'scripts.crust_review', ['--latest']),
    ('crust_budget', 'scripts.crust_budget', ['--latest']),
    ('elevation_check', 'scripts.elevation_check', ['--latest']),
    ('margin_width', 'scripts.margin_width', ['--latest']),
    ('interior_probe', 'scripts.interior_probe', ['--latest']),
)


def review_cmd(spec, sw, ctl, ref):
    """The review's command line, or None when it needs a control or an identity reference that was not given."""
    name, mod, args = spec
    if ('{REF}' in args and ref is None) or (ctl is None and ('{CTL}' in args or '{CTLNAME}' in args)):
        return None
    sub = {'{SW}': str(sw), '{CTL}': str(sw/ctl) if ctl else '', '{CTLNAME}': ctl or '', '{REF}': str(ref or ''),
           '{OUT}': str(sw/'reviews'/f'{name}.json')}
    argv = [sub.get(x, x) for x in args]
    head = [sys.executable, str(Path(mod))] if mod.endswith('.py') else [sys.executable, '-u', '-m', mod]
    return head+argv


def run_watched(cmd, out_file, tracker, detail, pred, pytest_progress=False):
    """Run a command with its output in out_file; the progress line follows its time (or pytest's [ NN%])."""
    t0 = time.time()
    with open(out_file, 'w', encoding='utf-8', errors='replace') as f:
        p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        while p.poll() is None:
            time.sleep(.5)
            frac = None
            if pytest_progress:
                try:
                    txt = Path(out_file).read_text(encoding='utf-8', errors='replace')
                    k = txt.rfind('%]')
                    if k > 0:
                        frac = int(txt[txt.rfind('[', 0, k)+1:k].strip())/100.
                except (OSError, ValueError):
                    pass
            spent = time.time()-t0
            if frac is None:
                frac = min(spent/pred, .99) if pred > 0 else 0.
            tracker.update(frac, detail, left=(spent*(1-frac)/frac if pytest_progress and frac > .1 else None))
    return p.returncode, time.time()-t0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--name', required=True, help='short sweep name (for the pytest printout and the timing record)')
    ap.add_argument('--grid', action='append', default=[])
    ap.add_argument('--seeds', default='1001,1002,1003,1004,1005,1006')
    ap.add_argument('--years', default='250')
    ap.add_argument('--rule', default='shelf')
    ap.add_argument('--save-every', default='5')
    ap.add_argument('--workers', default=None)
    ap.add_argument('--tests', nargs='*', default=None, help='test files (default: PLANET_TESTS)')
    ap.add_argument('--control', default=None, help="the control arm's folder name")
    ap.add_argument('--identity', default=None, help='folder the control must equal (tools/compare_runs.py)')
    ap.add_argument('--reviews', nargs='*', default=None, help='review names (default: all of REVIEWS that apply)')
    a = ap.parse_args(argv)
    from scripts import sweep_plates as spl
    cal = OUTPUTS/'planet-calibration'
    cal.mkdir(parents=True, exist_ok=True)
    tests = a.tests if a.tests is not None else [str(Path('tests')/f'{t}.py') for t in PLANET_TESTS]
    tests = [t for t in tests if Path(t).exists()]
    sweep_argv = ['--rule', a.rule, '--seeds', a.seeds, '--years', a.years, '--save-every', a.save_every]
    for gspec in a.grid:
        sweep_argv += ['--grid', gspec]
    if a.workers:
        sweep_argv += ['--workers', a.workers]
    workers = int(a.workers) if a.workers else max(1, min(8, (__import__('os').cpu_count() or 2)//2))
    n_settings = 1
    for gspec in a.grid:
        n_settings *= len(gspec.split('=', 1)[1].split(','))
    n_runs = n_settings*len(a.seeds.split(','))
    reviews = [r for r in REVIEWS if a.reviews is None or r[0] in a.reviews]
    per_run = pg.predict_run(workers, float(a.years))
    pred = {'tests': pg.predict('pytest') if tests else 0.,
            'sweep': pg.schedule_left([], n_runs, per_run, workers),
            'reviews': sum(pg.predict('review', r[0], per='runs')*n_runs for r in reviews)}
    tr = pg.Tracker([('tests', pred['tests']), ('sweep', pred['sweep']), ('reviews', pred['reviews'])])
    tr.log(f"{a.name}: {n_runs} runs ({n_settings} settings x {len(a.seeds.split(','))} seeds) on {workers} workers; "
           f"predicted ~{pg.fmt(sum(pred.values()))} (tests ~{pg.fmt(pred['tests'])}, runs ~{pg.fmt(pred['sweep'])}, "
           f"reviews ~{pg.fmt(pred['reviews'])})")
    # 1 tests
    tr.start_step()
    ptxt = cal/f'{a.name}-pytest.txt'
    if tests:
        rc, sec = run_watched([sys.executable, '-m', 'pytest', '-q', *tests], ptxt, tr, f'{len(tests)} test files',
                              pred['tests'], pytest_progress=True)
        last = [l for l in ptxt.read_text(encoding='utf-8', errors='replace').splitlines() if l.strip()][-1:]
        if rc != 0:
            tr.finish()
            print(ptxt.read_text(encoding='utf-8', errors='replace'))
            print(f'FAILED: the tests (printout above, also in {ptxt})')
            return 1
        pg.record('pytest', a.name, sec, files=len(tests))
        tr.log(f"tests passed: {last[0].strip() if last else ''} ({pg.fmt(sec)})")
    tr.end_step()
    # 2 sweep
    tr.start_step()
    sw = spl.main(sweep_argv, tracker=tr)
    tr.end_step()
    # 3 reviews
    tr.start_step()
    rdir = sw/'reviews'
    rdir.mkdir(exist_ok=True)
    if ptxt.exists():
        shutil.copyfile(ptxt, rdir/'pytest.txt')
    todo = [(r, review_cmd(r, sw, a.control, a.identity)) for r in reviews]
    todo = [(r, c) for r, c in todo if c is not None]
    failed = []
    preds = [pg.predict('review', r[0], per='runs')*n_runs for r, _ in todo]
    for k, ((r, cmd), p) in enumerate(zip(todo, preds)):
        base = sum(preds[:k])/max(sum(preds), 1e-9)
        share = p/max(sum(preds), 1e-9)

        class Sub:                                      # maps this review's own share onto the reviews step
            def update(self, frac=None, detail=None, left=None):
                tr.update(base+share*(frac or 0.), detail,
                          left=None if left is None else left+sum(preds[k+1:]))
        rc, sec = run_watched(cmd, rdir/f'{r[0]}.txt', Sub(), f'{r[0]} ({k+1} of {len(todo)})', p)
        pg.record('review', r[0], sec, runs=n_runs, sweep=sw.name)
        if rc != 0:
            failed.append(r[0])
            tr.log(f'review {r[0]} FAILED ({pg.fmt(sec)}): see {rdir/(r[0]+".txt")}')
        else:
            tr.log(f'review {k+1} of {len(todo)} done: {r[0]} ({pg.fmt(sec)})')
    tr.end_step()
    pg.record('total', a.name, tr.elapsed(), predicted=tr.predicted_at_start, runs=n_runs, workers=workers, sweep=sw.name)
    (rdir/'timings.txt').write_text('\n'.join(tr.summary())+'\n', encoding='utf-8')
    tr.finish(f"Done in {pg.fmt(tr.elapsed())} (predicted at the start ~{pg.fmt(tr.predicted_at_start)}). "
              f"Reviews in {rdir}" + (f"; FAILED: {', '.join(failed)}" if failed else ''))
    return 1 if failed else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except BaseException:
        sys.stdout.write('\n')                         # a traceback starts below the progress line
        raise
