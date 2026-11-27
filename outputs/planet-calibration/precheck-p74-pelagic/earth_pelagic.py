"""Pass 74: Earth's pelagic accumulation from GlobSed against floor age (Seton 2020) on n = 217, on open-ocean nodes
(label-free ocean, far from the continent and from Bird's trenches): the slope of bulk thickness with age."""
import sys
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid.paths import OUTPUTS
from planet_reviews.ocean_sediment_probe import earth_fields, distances
g = CubeSphere(217)
E = earth_fields(g)
dt_, dc_ = distances(g, E)
A = np.load(OUTPUTS/'planet-reference'/'earth-seafloor-age-n217.npz')
age, ds = A['age_myr'], A['dated_share']
sed = E['sed']
edges = np.arange(0, 181, 10.)
for cmin, tmin in ((200, 200), (500, 300), (800, 500), (1200, 500)):
    m = E['ocean'] & (ds > .5) & np.isfinite(age) & np.isfinite(sed) & (dc_ >= cmin) & (dt_ >= tmin)
    xs, med, mean, ar = [], [], [], []
    for a0, a1 in zip(edges[:-1], edges[1:]):
        k = m & (age >= a0) & (age < a1)
        if g.area[k].sum() < 1e12:   # < 1 Mkm2
            continue
        xs.append(np.average(age[k], weights=g.area[k])); med.append(np.median(sed[k]))
        mean.append(np.average(sed[k], weights=g.area[k])); ar.append(g.area[k].sum()/1e12)
    xs, med, mean, ar = map(np.array, (xs, med, mean, ar))
    w = ar
    def fit(y):
        X = np.c_[np.ones_like(xs), xs]; W = np.sqrt(w)[:, None]
        c = np.linalg.lstsq(X*W, y*W[:, 0], rcond=None)[0]; return c
    cm, cn = fit(med), fit(mean)
    print(f'continent >= {cmin} km, trench >= {tmin} km: area {g.area[m].sum()/1e12:.1f} Mkm2, mean age {np.average(age[m], weights=g.area[m]):.0f} Myr, '
          f'sediment mean {np.average(sed[m], weights=g.area[m]):.0f} m')
    print(f'   slope (area-weighted fit over bins) median {cm[1]:.2f} m/Myr (intercept {cm[0]:.0f} m), mean {cn[1]:.2f} m/Myr (intercept {cn[0]:.0f} m)')
    print('   bins (age: Mkm2 median/mean m):', ' | '.join(f'{x:.0f}: {a:.1f} {md:.0f}/{mn:.0f}' for x, a, md, mn in zip(xs, ar, med, mean)))

# The model's rate: grain (solid) thickness per Myr (sediment_compaction: sed is grain; the ocean's Athy column
# phi0 0.63, c 0.51 /km, PARAMS). Each node's GlobSed bulk h -> grain S = h - phi0/c (1 - exp(-c h)); the rate that
# gives Earth's far-ocean stock on its own floor ages: sum(a S) / sum(a age).
phi0, c = .63, .51e-3
print('grain rate (sum a S / sum a age), Athy phi0 0.63, c 0.51 /km:')
for cmin, tmin in ((200, 200), (500, 300), (800, 500), (1200, 500), (1600, 800)):
    m = E['ocean'] & (ds > .5) & np.isfinite(age) & np.isfinite(sed) & (dc_ >= cmin) & (dt_ >= tmin) & (age > 0)
    h = np.maximum(sed[m], 0); S = h-phi0/c*(1-np.exp(-c*h)); a = g.area[m]
    print(f'   continent >= {cmin}, trench >= {tmin} km: area {a.sum()/1e12:.1f} Mkm2, bulk {np.average(h, weights=a):.0f} m, '
          f'grain {np.average(S, weights=a):.0f} m, age {np.average(age[m], weights=a):.1f} Myr -> {(a*S).sum()/(a*age[m]).sum():.3f} m/Myr grain '
          f'({(a*h).sum()/(a*age[m]).sum():.2f} bulk)')
