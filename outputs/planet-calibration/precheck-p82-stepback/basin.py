"""Pass 82 step-back, probe 3: why the sea stands below the datum. The room below the datum (m of water per unit of
the whole surface, sum of area share x depth below the datum) by class, planet (six seeds) beside Earth (ETOPO; Seton
2020 ages; CRUST1.0 label-free continent, total >= 20 km). Classes: continent under the datum; ocean floor by age;
crust < 20 km without a floor age (Earth: undated in Seton's grid: margins, back-arcs; planet: continental-flag crust
thinned below 20 km). Per class the difference is split into area and depth (each at the other's mean, so the two
sum exactly). Then the ocean floor's depth by term on both (the model's law on Earth's age, CRUST1.0 crust, GlobSed
sediment; what is left is Earth's residual). Sea-level effect of a class: its excess room / the sea's own area share."""
import sys, numpy as np, common as C
from planet_tectonics.isostasy import gdh1_depth, bulk_thickness
arm = sys.argv[1] if len(sys.argv) > 1 else C.ARM
AG = [(0, 10), (10, 30), (30, 60), (60, 100), (100, 150), (150, 400)]
names = ['continent >= 20 km'] + [f'floor {a}-{b} Myr' for a, b in AG] + ['thin, no floor age']
e = C.earth(); p = e['p']; rc, rm, rw = p['rho_crust'], p['rho_mantle'], p['rho_water']
def classes(z, cont, oa):
    cl = [cont & (z < 0)] + [~cont & (z < 0) & (oa >= a) & (oa < b) for a, b in AG] + [~cont & (z < 0) & ~np.isfinite(oa)]
    return cl
def table(z, a, cl, terms=None):
    a = a/a.sum(); out = []
    for m in cl:
        ar = a[m].sum(); row = [ar, (C.wmean(-z[m], a[m]) if ar > 0 else 0.)]
        if terms is not None: row += [(C.wmean(t[m], a[m]) if ar > 0 else 0.) for t in terms]
        out.append(row)
    return np.array(out)
# Earth
ze = e['z']; conte = e['T'] >= 20e3
ds = np.load('outputs/planet-reference/earth-seafloor-age-n217.npz')['dated_share']   # the age grid holds 0 where undated
oae = np.where(conte | (ds < .5), np.nan, e['oa'])
h = np.where(np.isfinite(e['gsed']), e['gsed'], e['S']); phi0, c = p['phi0_ocean'], p['compaction_c_ocean_per_km']/1e3
solid = h-phi0/c*(1-np.exp(-c*h)); sed_e = h-(p['rho_grain']-rw)*solid/(rm-rw)
crust_e = (e['Tc']-p['ocean_thick_m'])*(rm-rc)/(rm-rw); gdh_e = -gdh1_depth(np.nan_to_num(oae))
res_e = ze-(gdh_e+crust_e+sed_e)
TE = table(ze, e['a'], classes(ze, conte, oae), [gdh_e, crust_e, sed_e, 0*ze, res_e, h, e['Tc']/1e3, np.nan_to_num(oae)])
print(f"Earth: room below the datum {np.sum(TE[:, 0]*TE[:, 1]):.0f} m (water 1.335e9 km3 = {1.335e18/e['a'].sum():.0f} m); surface below the datum {TE[:, 0].sum():.3f}")
for tag in ('t+0000.0', 't-0250.0'):
    TP = []; seas = []; asea = []
    for s in C.SEEDS:
        P = C.planet(s, tag, arm); z = P['zd']; contp = P['tot'] >= 20
        TP.append(table(z, P['a'], classes(z, contp, P['oa']), [P['gdh'], P['crust_o'], P['sedo'], P['dyn'], P['rest'], P['sh'], P['T']/1e3, P['oa']]))
        seas.append(P['sea']); asea.append(C.wmean((z < P['sea']).astype(float), P['a']))
    TPm = np.mean(TP, 0); seas = np.array(seas); A = float(np.mean(asea))
    roomP = np.array([np.sum(t[:, 0]*t[:, 1]) for t in TP]); roomE = np.sum(TE[:, 0]*TE[:, 1])
    print(f'\n=== {tag}: planet room below the datum {roomP.mean():.0f} m (sd {roomP.std(ddof=1):.0f}); excess over Earth {roomP.mean()-roomE:.0f} m; '
          f'sea {seas.mean():.0f} m; surface under the sea {A:.3f}; excess / that share = {-(roomP.mean()-roomE)/A:.0f} m')
    print('  class            | area share E / P | mean depth E / P (m) | room E / P (m) | excess: by area + by depth = total (sd seeds) | sea effect (m)')
    tot = 0
    for k, n in enumerate(names):
        aE, dE = TE[k, :2]; aP = np.array([t[k, 0] for t in TP]); dP = np.array([t[k, 1] for t in TP])
        byA = (aP-aE)*(dP+dE)/2; byD = (dP-dE)*(aP+aE)/2; ex = aP*dP-aE*dE; tot += ex.mean()
        print(f'  {n:18} | {aE:.3f} / {aP.mean():.3f} | {dE:5.0f} / {dP.mean():5.0f} | {aE*dE:5.0f} / {(aP*dP).mean():5.0f} | '
              f'{byA.mean():5.0f} + {byD.mean():5.0f} = {ex.mean():5.0f} ({ex.std(ddof=1):3.0f}) | {-ex.mean()/A:5.0f}')
    print(f'  sum of excess {tot:.0f} m')
    print('  ocean floor by term (area means, m): gdh | crust | sediment | dyn | rest (Earth: its residual) || sediment bulk m | crust km | age Myr')
    for k, n in enumerate(names[1:-1], 1):
        print(f'  {n:18} E: ' + ' | '.join(f'{v:6.0f}' for v in TE[k, 2:7]) + ' || ' + ' | '.join(f'{v:6.1f}' for v in TE[k, 7:]))
        print(f'  {"":18} P: ' + ' | '.join(f'{v:6.0f}' for v in TPm[k, 2:7]) + ' || ' + ' | '.join(f'{v:6.1f}' for v in TPm[k, 7:]))
    k = len(names)-1
    print(f'  {names[k]:18} E: ' + ' | '.join(f'{v:6.0f}' for v in TE[k, 2:7]) + ' || ' + ' | '.join(f'{v:6.1f}' for v in TE[k, 7:]))
    print(f'  {"":18} P: ' + ' | '.join(f'{v:6.0f}' for v in TPm[k, 2:7]) + ' || ' + ' | '.join(f'{v:6.1f}' for v in TPm[k, 7:]))
