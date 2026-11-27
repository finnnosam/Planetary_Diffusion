"""E3 step 4: E3's present-day field against E1's slab-proxy field (SMEAN2 fast part > 300 km) through the same kernels:
correlation per degree and over degrees 1-8, and the rms ratio."""
import json
import sys
import numpy as np
from planet_tectonics import mantle_flow as mf
from planet_reviews import mantle_flow_earth as mfe

d = np.load(sys.argv[1]); H3 = d['H']; lmax = int(d['lmax'])
out = {}
for name, key in (('E1 slab proxy (SMEAN2 fast > 300 km)', ('smean2', True, 300.)), ('E1 all SMEAN2 > 300 km', ('smean2', False, 300.))):
    l1, H1, _ = mfe.earth_flow(*key)
    assert l1 == lmax
    l, m = mf.alm_l(lmax), mf.alm_m(lmax)
    w = np.where(m == 0, 1., 2.)
    def corr(L):
        s = np.isin(l, L)
        num = (w[s]*(H3[s]*np.conj(H1[s])).real).sum()
        return float(num/np.sqrt((w[s]*np.abs(H3[s])**2).sum()*(w[s]*np.abs(H1[s])**2).sum()))
    P3, P1 = mf.spectrum(H3, lmax), mf.spectrum(H1, lmax)
    out[name] = {'corr_l1_8': corr(range(1, 9)), 'corr_l1_31': corr(range(1, 32)),
                 'corr_by_degree': {k: corr([k]) for k in (1, 2, 3, 4, 5, 6, 8, 12)},
                 'rms_E3_m': float(np.sqrt(P3[1:].sum())), 'rms_E1_m': float(np.sqrt(P1[1:].sum()))}
print(json.dumps(out, indent=1))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
