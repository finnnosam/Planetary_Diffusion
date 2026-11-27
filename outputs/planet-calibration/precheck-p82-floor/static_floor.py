"""Pass 82 (journal 107): the static size of ocean_depth_rule 'earth' / 'earth_all' on an arm's 0 Myr states (default:
today's model, crust / shorten of sweep-20261112-203456): every ocean-floor node moved by GDH1(age) minus the rule's
depth, continental columns re-floated, the sea re-solved for the water volume (as precheck-p82-stepback/counter.py);
no trajectory response. Six seeds: mean (sd between seeds). Run from the repo root: static_floor.py [arm]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'precheck-p82-stepback'))
import common as C
from planet_tectonics.isostasy import gdh1_depth, ocean_depth, sea_level_for_volume
arm = sys.argv[1] if len(sys.argv) > 1 else 'outputs/planet-calibration/sweep-20261112-203456/start_lithosphere-crust_belt_lithosphere-shorten'
def resolve(P, zd_new):
    p = P['p']; rm, rw = p['rho_mantle'], p['rho_water']; V = p['water_volume_km3']*1e9
    a, cont, s0 = P['a'], P['cont'], P['sea']; f = rm/(rm-rw)
    dry = np.where(cont & (P['zd'] < s0), s0+(P['zd']-s0)/f, P['zd'])+(zd_new-P['zd']); s = s0
    for _ in range(30):
        z = np.where(cont & (dry < s), s+(dry-s)*f, dry); s1 = sea_level_for_volume(z, a, V)
        if abs(s1-s) < .01: break
        s = s1
    return z, s1
R = {k: [] for k in ('gdh1', 'earth', 'earth_all')}
for sd in C.SEEDS:
    P = C.planet(sd, 't+0000.0', arm); a = P['a']/P['a'].sum(); oc = ~P['cont']
    for k in R:
        lift = np.where(oc, gdh1_depth(np.nan_to_num(P['oa']))-ocean_depth(np.nan_to_num(P['oa']), {**P['p'], 'ocean_depth_rule': k}), 0.)
        z, s = resolve(P, P['zd']+lift); l = z > s
        R[k].append([s, a[l].sum(), *C.wq(z[l]-s, a[l], [.25, .5, .9]), np.sum(a*np.maximum(-z, 0)), np.sum(a[oc]*lift[oc])/a[oc].sum()])
base = np.array(R['gdh1'])
print(arm); print('rule | sea on the datum | land share | land p25 / p50 / p90 | room below the datum | mean lift of the floor || change of the sea (sd) | of p50 (sd) | of p90')
for k in R:
    r = np.array(R[k]); m = r.mean(0); d = r-base
    print(f'  {k:9} | {m[0]:5.0f} | {m[1]:.3f} | {m[2]:.0f} / {m[3]:.0f} / {m[4]:.0f} | {m[5]:.0f} | {m[6]:.0f} || {d[:, 0].mean():+.0f} ({d[:, 0].std(ddof=1):.0f}) | {d[:, 3].mean():+.0f} ({d[:, 3].std(ddof=1):.0f}) | {d[:, 4].mean():+.0f}')
