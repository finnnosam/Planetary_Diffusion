"""Pass 82 step-back, probe 12: the sea against the continents' own elevation distribution, both ways round.
(1) Earth's surface (ETOPO) given the planet's sea level on the datum: land p50 of Earth with the sea at -487 m.
(2) The planet with its sea at the datum (hyps.py). Their mean is each factor's share whatever the order.
(3) Within the continents (crust + sediment >= 30 km): every column shifted by Earth's minus the planet's median on the
datum (level only), or mapped rank for rank onto Earth's elevation distribution of the same class (level and spread);
the sea as run, re-solved for the water volume (plain volume balance), and at the datum."""
import numpy as np, common as C
from planet_tectonics.isostasy import sea_level_for_volume
e = C.earth(); ae = e['a']; ze = e['z']; V = e['p']['water_volume_km3']*1e9
de = (e['T'] >= 30e3) & ~e['ice']; QQ = np.linspace(0, 1, 2001); qe = C.wq(ze[de], ae[de], QQ)
def land(z, a, s):
    l = z > s; return [a[l].sum()/a.sum(), *C.wq(z[l]-s, a[l], [.25, .5, .9])]
R = {k: [] for k in ('as run', 'Earth, the planet\'s sea', 'planet, sea at the datum', 'level: sea as run', 'level: sea re-solved', 'level: sea at the datum',
                     'level+spread: sea as run', 'level+spread: sea re-solved', 'level+spread: sea at the datum')}
seas = {k: [] for k in R}
for s in C.SEEDS:
    P = C.planet(s); a = P['a']; z = P['zd']; sea = P['sea']; d = (P['tot'] >= 30) & P['cont']
    R['as run'].append(land(z, a, sea)); seas['as run'].append(sea)
    R["Earth, the planet's sea"].append(land(ze, ae, sea)); seas["Earth, the planet's sea"].append(sea)
    R['planet, sea at the datum'].append(land(z, a, 0.)); seas['planet, sea at the datum'].append(0.)
    z1 = z.copy(); z1[d] += C.wq(ze[de], ae[de], .5)-C.wq(z[d], a[d], .5)
    o = np.argsort(z[d]); cw = np.cumsum(a[d][o]); r = np.empty(d.sum()); r[o] = (cw-.5*a[d][o])/cw[-1]
    z2 = z.copy(); z2[d] = np.interp(r, QQ, qe)
    for nm, zz in (('level', z1), ('level+spread', z2)):
        s1 = sea_level_for_volume(zz, a, V)
        for lab, sv in (('sea as run', sea), ('sea re-solved', s1), ('sea at the datum', 0.)):
            R[f'{nm}: {lab}'].append(land(zz, a, sv)); seas[f'{nm}: {lab}'].append(sv)
le = ze > 0; E = [ae[le].sum()/ae.sum(), *C.wq(ze[le], ae[le], [.25, .5, .9])]
print(f'Earth: land {E[0]:.3f}, p25 / p50 / p90 {E[1]:.0f} / {E[2]:.0f} / {E[3]:.0f}')
print('scenario | sea on the datum | land share | land p25 / p50 / p90 | p50 sd between seeds')
for k in R:
    r = np.array(R[k]); m = r.mean(0); print(f'  {k:32} | {np.mean(seas[k]):5.0f} | {m[0]:.3f} | {m[1]:.0f} / {m[2]:.0f} / {m[3]:.0f} | {r[:, 2].std(ddof=1):.0f}')
A = np.mean([r[2] for r in R['as run']]); X = np.mean([r[2] for r in R["Earth, the planet's sea"]]); D = np.mean([r[2] for r in R['planet, sea at the datum']]); E50 = E[2]
print(f'\nexcess {A-E50:.0f} m. The sea: {A-D:.0f} (moved first) or {X-E50:.0f} (moved last), mean {(A-D+X-E50)/2:.0f} = {100*(A-D+X-E50)/2/(A-E50):.0f}%. '
      f'The continents on the datum: {A-X:.0f} (first) or {D-E50:.0f} (last), mean {(A-X+D-E50)/2:.0f} = {100*(A-X+D-E50)/2/(A-E50):.0f}%.')
