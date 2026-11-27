"""Pass 82 step-back, probe 13: three factors, every combination: S the sea (as run / at the datum), L the continents'
level (crust + sediment >= 30 km shifted to Earth's median on the datum), W their spread (mapped rank for rank onto
Earth's distribution about the median; W without L keeps the planet's median). Land p50 at each of the eight corners
and each factor's mean contribution over the six orders (Shapley). Usage: three_factor.py [arm] Then W split: only the columns below the
continent's median mapped (low side), or only those above (high side). No sea re-solve (S is the factor)."""
import itertools, sys, numpy as np, common as C
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM
e = C.earth(); ae = e['a']; ze = e['z']; de = (e['T'] >= 30e3) & ~e['ice']; QQ = np.linspace(0, 1, 2001); qe = C.wq(ze[de], ae[de], QQ); me = C.wq(ze[de], ae[de], .5)
def p50(z, a, s): l = z > s; return C.wq(z[l]-s, a[l], .5)
F = {}; side = {k: [] for k in ('low side, sea as run', 'high side, sea as run', 'low side, sea at the datum', 'high side, sea at the datum')}
for s in C.SEEDS:
    P = C.planet(s, 't+0000.0', arm); a = P['a']; z = P['zd']; sea = P['sea']; d = (P['tot'] >= 30) & P['cont']; mp = C.wq(z[d], a[d], .5)
    o = np.argsort(z[d]); cw = np.cumsum(a[d][o]); r = np.empty(d.sum()); r[o] = (cw-.5*a[d][o])/cw[-1]; zmap = np.interp(r, QQ, qe)
    for S, L, W in itertools.product((0, 1), repeat=3):
        zz = z.copy()
        if W: zz[d] = zmap-(0 if L else me-mp)
        elif L: zz[d] += me-mp
        F.setdefault((S, L, W), []).append(p50(zz, a, 0. if S else sea))
    for nm, sel in (('low side', r < .5), ('high side', r >= .5)):
        zz = z.copy(); zz[d] += me-mp; t = zz[d]; t[sel] = zmap[sel]; zz[d] = t
        side[f'{nm}, sea as run'].append(p50(zz, a, sea)); side[f'{nm}, sea at the datum'].append(p50(zz, a, 0.))
f = {k: float(np.mean(v)) for k, v in F.items()}
print(arm)
print('land p50 (m) at the corners (S sea at the datum, L level, W spread): ' + ', '.join(f'{"".join(n for n, on in zip("SLW", k) if on) or "none"} {v:.0f}' for k, v in sorted(f.items())))
sh = dict(S=0., L=0., W=0.)
for order in itertools.permutations('SLW'):
    on = dict(S=0, L=0, W=0)
    for x in order:
        before = f[(on['S'], on['L'], on['W'])]; on[x] = 1; sh[x] += (before-f[(on['S'], on['L'], on['W'])])/6
tot = f[(0, 0, 0)]-f[(1, 1, 1)]
print(f'excess closed by the three together {tot:.0f} m (Earth\'s p50 380; all three give {f[(1, 1, 1)]:.0f}). Mean contribution over the orders: sea {sh["S"]:.0f} ({100*sh["S"]/tot:.0f}%), '
      f'level {sh["L"]:.0f} ({100*sh["L"]/tot:.0f}%), spread {sh["W"]:.0f} ({100*sh["W"]/tot:.0f}%)')
rng = lambda x: f'{min(x):.0f} ... {max(x):.0f}'
alone = {'S': [f[(0, 0, 0)]-f[(1, 0, 0)], f[(0, 1, 0)]-f[(1, 1, 0)], f[(0, 0, 1)]-f[(1, 0, 1)], f[(0, 1, 1)]-f[(1, 1, 1)]],
         'L': [f[(0, 0, 0)]-f[(0, 1, 0)], f[(1, 0, 0)]-f[(1, 1, 0)], f[(0, 0, 1)]-f[(0, 1, 1)], f[(1, 0, 1)]-f[(1, 1, 1)]],
         'W': [f[(0, 0, 0)]-f[(0, 0, 1)], f[(1, 0, 0)]-f[(1, 0, 1)], f[(0, 1, 0)]-f[(0, 1, 1)], f[(1, 1, 0)]-f[(1, 1, 1)]]}
print('range of each factor\'s effect over what else is already changed: sea ' + rng(alone['S']) + ', level ' + rng(alone['L']) + ', spread ' + rng(alone['W']))
print('spread, one side at a time (level applied): ' + '; '.join(f'{k} {np.mean(v):.0f}' for k, v in side.items()) + f'   (level only: sea as run {f[(0, 1, 0)]:.0f}, at the datum {f[(1, 1, 0)]:.0f}; both sides: {f[(0, 1, 1)]:.0f}, {f[(1, 1, 1)]:.0f})')
