import sys, json, runpy
from planet_tectonics import events as E
from planet_tectonics.breakup import production_per_mkm2
log = []
orig = E.EventsMixin._split_margin_km
def wrapped(self, own, lab):
    L = orig(self, own, lab)
    if self.p['rift_line'] in ('force', 'chord'):
        A = self.g.area; Ac = A[own & self.world['cont']].sum()/1e12
        part = (lab == 1) & own
        log.append((float(self.t), int(self.owner[own][0]), getattr(self, '_line_kind', None), round(L), round(Ac, 1),
                    round(A[part].sum()/1e12, 1), production_per_mkm2(self.t-self.t_start)))
    return L
E.EventsMixin._split_margin_km = wrapped
sys.argv = ['simulate_plates', '--seeds', '1007', '--years', sys.argv[2], '--save-every', '1000', '--save-last', '0',
            '--out', sys.argv[3], '--param', 'threads=1', '--param', f'rift_line={sys.argv[1]}']
try:
    runpy.run_module('scripts.simulate_plates', run_name='__main__')
finally:
    json.dump(log, open(sys.argv[-1].split('=')[0] if False else '/home/claude/p73/Llog_%s.json' % sys.argv[-1].split('=')[1], 'w'))
