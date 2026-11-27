"""Pass 70: Earth's ridge-near continental edge (K1's instrument) split by Bird's class: CRB (continental rift, the
planet's stalled-rift analogue) and OSR (ocean ridge), with distance < 100 km."""
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_grid import fields
from planet_grid.paths import DATA
from planet_reviews.crust_review import label_free_continent
from planet_reviews.margin_width import earth_nodes
from planet_earth.bird import read_dbf, xyz_from_lonlat
from planet_tectonics.graphops import nearest_label
g = CubeSphere(217); nb = g.nbr; ok = nb >= 0; sp = float(np.sqrt(g.area.mean())); k = lambda m: m.sum()*sp/1e6
T, _ = earth_nodes(g); cont = label_free_continent(T)
s = read_dbf(DATA/'tectonicplates-master'/'PB2002_steps.dbf')
x = xyz_from_lonlat(s['STARTLONG'], s['STARTLAT'])+xyz_from_lonlat(s['FINALLONG'], s['FINALLAT'])
x /= np.linalg.norm(x, axis=1, keepdims=True)
bt = np.full(g.size, '', dtype='<U3'); i = fields.nearest_nodes(g, x); bt[i] = [c.strip() for c in s['STEPCLASS']]
vd = np.zeros(g.size); vd[i] = s['VELOCITYDI']
edge = cont & np.any(np.where(ok, ~cont[np.where(ok, nb, 0)], False), 1)
bl, bd = nearest_label(np.where(bt != '', np.arange(g.size), -1).astype(np.int64), nb, g.dist)
J = np.maximum(bl, 0); near = edge & (bd <= 300e3)
for c in ('CRB', 'OSR'):
    m = near & (bt[J] == c)
    print(f"Earth ridge-near via {c}: {k(m):.1f}k km (d < 100 km {k(m & (bd < 1e5)):.1f}k); opening p50 {np.median(vd[J][m]):.0f} mm/yr")
