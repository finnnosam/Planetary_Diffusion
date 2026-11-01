"""The save cadence (forty-eighth pass): the run records crust_budget's trench arrivals in its stats series, equal to
the arrivals computed from the saved snapshots, so sweeps can save every 5 Myr; compare_runs compares what both runs
have and lists the rest apart."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CRUST1 = ROOT/'data'/'crust1'


def _compare_runs():
    spec = importlib.util.spec_from_file_location('compare_runs', ROOT/'tools'/'compare_runs.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


@pytest.mark.skipif(not CRUST1.is_dir(), reason='CRUST1.0 not present')
def test_series_arrivals_match_snapshots_and_cadence_keeps_identity(tmp_path):
    from scripts import simulate_plates
    from planet_grid.sphere_graph import CubeSphere
    from planet_reviews import crust_budget as cb
    from planet_reviews.crust_review import run_params
    base = ['--n', '96', '--seeds', '7', '--years', '12', '--save-last', '6', '--param', 'threads=1']
    simulate_plates.main(base+['--save-every', '1', '--out', str(tmp_path/'a')])
    simulate_plates.main(base+['--save-every', '5', '--out', str(tmp_path/'b')])
    a, b = tmp_path/'a'/'seed7', tmp_path/'b'/'seed7'
    stats = json.loads((a/'stats.json').read_text())
    assert all('trench_arrivals' in r for r in stats['series'])
    arr = cb.series_arrivals(stats, -6, 0)
    assert arr is not None and len(arr) == 7
    g, p = CubeSphere(96), run_params(a)
    for t, x in zip(range(-6, 1), arr):
        y = cb.trench_arrivals(np.load(a/'snapshots'/f't{float(t):+07.1f}.npz'), g, p)
        assert x == pytest.approx(y, rel=0, abs=0)
    # a run before the pass (no series arrivals) falls back to the snapshots
    assert cb.series_arrivals({'series': [{'time_myr': float(t)} for t in range(-6, 1)]}, -6, 0) is None
    na = len(list((a/'snapshots').glob('t*.npz'))); nb = len(list((b/'snapshots').glob('t*.npz')))
    assert nb < na
    cr = _compare_runs()
    assert cr.compare(a, b)                    # the simulation does not depend on what it saves
    assert cr.compare(b, a)
