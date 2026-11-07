"""Progress display and timing record for sweeps (sixty-eighth pass, journal 90; the user asked on 7 Nov for elapsed
against predicted time, overall and step percentages, plain 'run k done' lines and one progress bar that updates in
place, and for every step's time to be recorded so later predictions stay accurate).

  Tracker   one line at the bottom of the console, redrawn in place:
              [#########-----------]  45% | 0:27:10 of ~1:00:20 | 2/3 sweep: 9 of 24 runs done, 8 running [41%]
            with event lines printed above it ('run 9 of 24 done: settings 1 of 4 (...), seed 1003 (3 of 6), 18.2
            min'). When the output is not a console (redirected to a file) it prints a plain line every 5% instead.
  record()  appends one JSON line per finished step to outputs/planet-calibration/timings.jsonl (kind, name, seconds
            and what it depends on: runs, workers, years).
  predict() the median of the last matching records, or DEFAULTS (measured on sweep-20261107-094624: 24 runs of
            250 Myr on 8 workers, 61 min; run median 19.6 min, pytest 116 s, the reviews 4.5 min).
Output only: nothing here touches a run's state.
"""
import json
import math
import os
import shutil
import statistics
import sys
import time

from planet_grid.paths import OUTPUTS

TIMINGS = 'timings.jsonl'
# seconds; run: one 250 Myr n = 217 run while 8 run at once; review: per run reviewed (the reviews read every run)
DEFAULTS = {('pytest', None): 120., ('run', None): 21*60.,
            ('review', 'identity'): .6, ('review', 'max_elevation'): .2, ('review', 'rift_production'): .5,
            ('review', 'passive_survival'): .5, ('review', 'coast_budget'): 1.4, ('review', 'coast_profile'): 1.2,
            ('review', 'plate_pieces_review'): 4.1, ('review', 'crust_review'): .3, ('review', 'crust_budget'): 1.3,
            ('review', 'elevation_check'): .5, ('review', 'margin_width'): .6, ('review', 'interior_probe'): .6}
RECENT = 8


def timings_path():
    return OUTPUTS/'planet-calibration'/TIMINGS


def fmt(seconds):
    """h:mm:ss"""
    s = int(round(max(seconds, 0.)))
    return f'{s//3600}:{s % 3600//60:02d}:{s % 60:02d}'


def record(kind, name=None, seconds=0., path=None, **meta):
    """Append one timing record (output only)."""
    p = path or timings_path()
    try:
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, 'a', encoding='utf-8') as f:
            f.write(json.dumps({'when': time.strftime('%Y-%m-%d %H:%M:%S'), 'kind': kind, 'name': name,
                                'seconds': float(seconds), **meta}, default=str)+'\n')
    except OSError:
        pass                                             # a timing record never stops a sweep


def read_records(path=None):
    p = path or timings_path()
    out = []
    try:
        for line in open(p, encoding='utf-8'):
            try:
                out.append(json.loads(line))
            except ValueError:
                continue
    except OSError:
        pass
    return out


def predict(kind, name=None, per=None, path=None, **match):
    """Seconds for one step: the median of the last RECENT records of this kind (and name) whose `match` fields are
    equal; with `per` given, records are divided by their own value of that field (seconds per run) and the result is
    per unit. Falls back to DEFAULTS."""
    rs = [r for r in read_records(path) if r.get('kind') == kind and (name is None or r.get('name') == name)
          and all(r.get(k) == v for k, v in match.items())]
    if per is not None:
        rs = [r for r in rs if r.get(per)]
    if rs:
        vals = [r['seconds']/(r[per] if per else 1.) for r in rs[-RECENT:]]
        return float(statistics.median(vals))
    return DEFAULTS.get((kind, name), DEFAULTS.get((kind, None), 60.))


def predict_run(workers, years, path=None):
    """Seconds for one run: recent runs with the same workers and length; else recent runs of any setup scaled by
    length; else the default scaled by length."""
    rs = [r for r in read_records(path) if r.get('kind') == 'run' and r.get('years')]
    same = [r['seconds'] for r in rs if r.get('workers') == workers and r.get('years') == years]
    if same:
        return float(statistics.median(same[-RECENT*4:]))
    if rs:
        return float(statistics.median([r['seconds']*years/r['years']
                                        for r in rs[-RECENT*4:]]))
    return DEFAULTS[('run', None)]*years/250.


def schedule_left(running_left, queued, per_run, workers):
    """Seconds until a pool of `workers` finishes: the running runs' remaining seconds, then `queued` runs of
    `per_run` seconds each, every one going to the worker free first."""
    free = sorted(list(running_left)+[0.]*max(0, workers-len(running_left)))[:max(workers, len(running_left))]
    if not free:
        free = [0.]
    for _ in range(queued):
        free.sort()
        free[0] += per_run
    return max(free)


class Tracker:
    """Steps with predicted durations; one progress line for the whole, the current step's share and detail."""

    def __init__(self, steps, stream=None, live=None, width=None):
        """steps: [(name, predicted seconds)] in order."""
        self.steps = [[n, float(p), None] for n, p in steps]   # name, predicted, actual
        self.stream = stream or sys.stdout
        self.live = (self.stream.isatty() if hasattr(self.stream, 'isatty') else False) if live is None else live
        self.width = width
        self.t0 = time.time()
        self.k = -1
        self.step_t0 = self.t0
        self.step_frac, self.step_left, self.detail = 0., None, ''
        self.predicted_at_start = sum(p for _, p, _ in self.steps)
        self.shown = False
        self.last_plain = -1
        self.last_draw = 0.

    # ---- steps
    def start_step(self, name=None, predicted=None):
        self.k += 1
        if name is not None:
            self.steps[self.k][0] = name
        if predicted is not None:
            self.steps[self.k][1] = float(predicted)
        self.step_t0 = time.time()
        self.step_frac, self.step_left, self.detail = 0., None, ''
        self.draw(force=True)

    def end_step(self):
        self.steps[self.k][2] = time.time()-self.step_t0
        self.step_frac, self.step_left = 1., 0.
        self.draw(force=True)
        return self.steps[self.k][2]

    def update(self, frac=None, detail=None, left=None):
        """The current step's fraction done (0-1), a short detail, and (optional) its remaining seconds if known
        better than from its prediction."""
        if frac is not None:
            self.step_frac = min(max(frac, 0.), 1.)
        if detail is not None:
            self.detail = detail
        self.step_left = left
        self.draw()

    # ---- estimates
    def elapsed(self):
        return time.time()-self.t0

    def current_left(self):
        if self.k < 0 or self.k >= len(self.steps):
            return 0.
        if self.step_left is not None:
            return self.step_left
        spent = time.time()-self.step_t0
        pred = self.steps[self.k][1]
        if self.step_frac > .1:                          # far enough in: the step's own pace
            return spent*(1-self.step_frac)/self.step_frac
        return max(pred-spent, pred*(1-self.step_frac), 0.)

    def total(self):
        """Predicted whole: elapsed + this step's remainder + the steps not started."""
        later = sum(p for _, p, _ in self.steps[self.k+1:])
        return self.elapsed()+self.current_left()+later

    def overall(self):
        tot = self.total()
        return self.elapsed()/tot if tot > 0 else 1.

    # ---- drawing
    def line(self):
        w = self.width or shutil.get_terminal_size((100, 20)).columns
        f = self.overall()
        nb = 20
        bar = '#'*int(round(f*nb))+'-'*(nb-int(round(f*nb)))
        name = self.steps[self.k][0] if 0 <= self.k < len(self.steps) else ''
        s = (f"[{bar}] {f:4.0%} | {fmt(self.elapsed())} of ~{fmt(self.total())} | "
             f"{self.k+1}/{len(self.steps)} {name}: {self.detail} [{self.step_frac:.0%}]")
        return s[:max(w-1, 20)]

    def draw(self, force=False):
        if self.k < 0:                                   # nothing to show before the first step
            return
        now = time.time()
        if not force and now-self.last_draw < .5:
            return
        self.last_draw = now
        if self.live:
            s = self.line()
            self.stream.write('\r'+s+' '*max(0, (self.width or shutil.get_terminal_size((100, 20)).columns)-1-len(s)))
            self.stream.flush()
            self.shown = True
        else:
            pct = int(self.overall()*100)//5*5
            if pct != self.last_plain or force:
                self.last_plain = pct
                self.stream.write(self.line()+'\n'); self.stream.flush()

    def log(self, msg):
        """An event line above the progress line."""
        if self.live and self.shown:
            w = self.width or shutil.get_terminal_size((100, 20)).columns
            self.stream.write('\r'+' '*(w-1)+'\r')
        self.stream.write(msg+'\n')
        self.stream.flush()
        self.draw(force=True)

    def finish(self, msg=None):
        if self.live and self.shown:
            self.stream.write('\n')
        if msg:
            self.stream.write(msg+'\n')
        self.stream.flush()

    def summary(self):
        """Lines: each step's predicted and actual time, and the whole."""
        rows = [f"{'step':32s} {'predicted':>10s} {'actual':>10s}"]
        for n, p, a in self.steps:
            rows.append(f"{n:32s} {fmt(p):>10s} {fmt(a) if a is not None else '-':>10s}")
        rows.append(f"{'total':32s} {fmt(self.predicted_at_start):>10s} {fmt(self.elapsed()):>10s}")
        return rows
