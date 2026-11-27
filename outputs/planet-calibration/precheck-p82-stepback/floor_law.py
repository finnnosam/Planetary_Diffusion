"""Pass 82 step-back, probe 5: does the model's ocean-floor law hold on Earth? ETOPO on dated floor (Seton 2020, node
>= 50% dated, CRUST1.0 total < 20 km) minus GDH1(age) + the model's sediment term on GlobSed + its crust term on
CRUST1.0: Earth's residual by age class: area mean, quantiles, and by CRUST1.0 type (A1 normal oceanic, A0 young,
B- melt-affected / plateaus, others), and by distance from the label-free continent. Then the sea budget's area part
split into the deep floor's total area and its age mix."""
import numpy as np, common as C
from planet_tectonics.isostasy import gdh1_depth
from planet_landscape import graph_landscape as gl
e = C.earth(); p = e['p']; g = e['g']; rc, rm, rw = p['rho_crust'], p['rho_mantle'], p['rho_water']; a = e['a']/e['a'].sum()
ds = np.load('outputs/planet-reference/earth-seafloor-age-n217.npz')['dated_share']
fl = (e['T'] < 20e3) & (ds >= .5) & (e['z'] < 0); oa = e['oa']
h = e['gsed']; phi0, c = p['phi0_ocean'], p['compaction_c_ocean_per_km']/1e3
sed = h-(p['rho_grain']-rw)*(h-phi0/c*(1-np.exp(-c*h)))/(rm-rw); cr = (e['Tc']-p['ocean_thick_m'])*(rm-rc)/(rm-rw)
res = e['z']-(-gdh1_depth(oa)+cr+sed)
dc = gl.distance_to(e['T'] >= 20e3, g.nbr, g.dist)/1e3
AG = [(0, 10), (10, 30), (30, 60), (60, 100), (100, 150), (150, 400)]
print('Earth residual (ETOPO - GDH1 - sediment - crust terms), m. age class | area | mean | p10 / p25 / p50 / p75 / p90 | share > +500 m | mean by type A1 / A0 / B- / other (area share) | mean by distance from the continent < 300 / 300-1000 / > 1000 km')
for lo, hi in AG:
    m = fl & (oa >= lo) & (oa < hi); q = C.wq(res[m], a[m], [.1, .25, .5, .75, .9])
    ty = []
    for name, sel in (('A1', e['tn'] == 'A1'), ('A0', e['tn'] == 'A0'), ('B-', e['tn'] == 'B-'), ('other', ~np.isin(e['tn'], ['A1', 'A0', 'B-']))):
        k = m & sel; ty.append(f'{C.wmean(res[k], a[k]):.0f} ({a[k].sum()/a[m].sum():.2f})' if k.sum() > 20 else '-')
    dd = [f'{C.wmean(res[m & s], a[m & s]):.0f}' if (m & s).sum() > 20 else '-' for s in (dc < 300, (dc >= 300) & (dc < 1000), dc >= 1000)]
    print(f'  {lo:3d}-{hi:3d} | {a[m].sum():.3f} | {C.wmean(res[m], a[m]):5.0f} | ' + ' / '.join(f'{v:.0f}' for v in q) + f' | {a[m & (res > 500)].sum()/a[m].sum():.2f} | ' + ' / '.join(ty) + ' | ' + ' / '.join(dd))
m = fl & (oa >= 10); far = m & (dc >= 1000) & (e['tn'] == 'A1')
print(f'all dated floor >= 10 Myr: mean {C.wmean(res[m], a[m]):.0f}, median {C.wq(res[m], a[m], .5):.0f}; normal crust (A1) >= 1,000 km from the continent: mean {C.wmean(res[far], a[far]):.0f}, median {C.wq(res[far], a[far], .5):.0f} (area {a[far].sum():.3f})')
# the planet's own floor residual for comparison (dyn + rest), and the area / age-mix split
EA = np.array([a[fl & (oa >= lo) & (oa < hi)].sum() for lo, hi in AG]); ED = np.array([C.wmean(-e['z'][fl & (oa >= lo) & (oa < hi)], a[fl & (oa >= lo) & (oa < hi)]) for lo, hi in AG])
out = []
for s in C.SEEDS:
    P = C.planet(s); ap = P['a']/P['a'].sum(); z = P['zd']; f = (P['tot'] < 20) & (z < 0) & np.isfinite(P['oa'])
    PA = np.array([ap[f & (P['oa'] >= lo) & (P['oa'] < hi)].sum() for lo, hi in AG]); PD = np.array([C.wmean(-z[f & (P['oa'] >= lo) & (P['oa'] < hi)], ap[f & (P['oa'] >= lo) & (P['oa'] < hi)]) for lo, hi in AG])
    Dm = (PD+ED)/2; tot = (PA.sum()-EA.sum())*np.sum(Dm*(PA/PA.sum()+EA/EA.sum())/2); mix = np.sum((PA/PA.sum()-EA/EA.sum())*Dm)*(PA.sum()+EA.sum())/2
    out.append([PA.sum(), tot, mix, np.sum((PA-EA)*Dm), np.sum((PD-ED)*(PA+EA)/2), C.wmean((z < P['sea']).astype(float), ap), np.sum(PA*P and 0 or 0) if False else 0])
o = np.array(out)[:, :6]; m_, s_ = o.mean(0), o.std(0, ddof=1); A = m_[5]
print(f'\ndated deep floor, share of the surface: Earth {EA.sum():.3f}, planet {m_[0]:.3f} (sd {s_[0]:.3f})')
print(f'room below the datum (m): more floor area {m_[1]:.0f} (sd {s_[1]:.0f}) + younger age mix {m_[2]:.0f} ({s_[2]:.0f}) = area part {m_[3]:.0f}; depth at a given age {m_[4]:.0f} ({s_[4]:.0f})')
print(f'as sea level (/ {A:.3f}): floor area {-m_[1]/A:.0f}, age mix {-m_[2]/A:.0f}, depth at a given age {-m_[4]/A:.0f} m')
# the depth part by term: per age class (planet term mean - Earth term mean) x mean area share, summed; as sea level
T = []
for s in C.SEEDS:
    P = C.planet(s); ap = P['a']/P['a'].sum(); z = P['zd']; f = (P['tot'] < 20) & (z < 0) & np.isfinite(P['oa']); row = np.zeros(5)
    for k, (lo, hi) in enumerate(AG):
        mp = f & (P['oa'] >= lo) & (P['oa'] < hi); me = fl & (oa >= lo) & (oa < hi); w = (ap[mp].sum()+a[me].sum())/2
        pe = [C.wmean(t[mp], ap[mp]) for t in (P['gdh'], P['sedo'], P['crust_o'], P['dyn']+P['rest'])]
        ee = [C.wmean(t[me], a[me]) for t in (-gdh1_depth(oa), sed, cr, res)]
        row += w*np.array([pe[0]-ee[0], pe[1]-ee[1], pe[2]-ee[2], pe[3], -ee[3]])
    T.append(row)
T = np.array(T); tm, ts = T.mean(0)/A, T.std(0, ddof=1)/A
print('depth at a given age, as sea level (m; sd seeds): ' + ', '.join(f'{n} {v:+.0f} ({s:.0f})' for n, v, s in zip(
    ['age within the class (GDH1)', 'sediment on the floor', 'ocean crust thickness', "the planet's own dyn + flexure", "Earth's floor above the law (its residual)"], tm, ts)) + f'; sum {tm.sum():+.0f}')
