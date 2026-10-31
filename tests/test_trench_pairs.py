"""planet_reviews.trench_pairs (forty-fourth pass): plate lineage, the relative-motion ratio and the coherence measure."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_reviews import trench_pairs as tp

PAR = {12: (8, 'rift', -230.), 23: (12, 'rift', -119.), 19: (3, 'rift', -162.),
       15: (6, 'subduction initiation', -167.), 16: (15, 'subduction initiation', -165.),
       13: (0, 'rift', -199.), 14: (0, 'rift', -183.), 25: (5, 'detach', -94.)}


def test_lineage_reads_split_events_only():
    ev = [{'t': -230., 'kind': 'rift', 'plate': 8, 'new_plate': 12, 'area_mkm2': 27.},
          {'t': -210., 'kind': 'merge', 'plate': 4, 'merged': 2, 'area_mkm2': 4.6},
          {'t': -94., 'kind': 'detach', 'plate': 5, 'new_plate': 25, 'area_mkm2': 2.9}]
    assert tp.lineage(ev) == {12: (8, 'rift', -230.), 25: (5, 'detach', -94.)}


def test_pair_relation():
    assert tp.pair_relation(19, 23, PAR) == ('different start plates', -250.)
    assert tp.pair_relation(8, 23, PAR) == ('rift', -230.)          # 23 descends from 8 through 12
    assert tp.pair_relation(23, 12, PAR) == ('rift', -119.)
    assert tp.pair_relation(6, 16, PAR) == ('subduction initiation', -167.)
    assert tp.pair_relation(15, 16, PAR) == ('subduction initiation', -165.)
    assert tp.pair_relation(13, 14, PAR) == ('rift', -199.)          # siblings: separate since the first split
    assert tp.pair_relation(25, 16, PAR) == ('different start plates', -250.)


def test_rho():
    v = np.array([[1., 0, 0], [1., 0, 0], [0, 2., 0]])
    w = np.array([[1., 0, 0], [-1., 0, 0], [0, 1., 0]])
    assert np.allclose(tp.rho(v, w), [0., 2., .5])


def test_coherence_signs():
    g = CubeSphere(16)
    pid = (g.xyz[:, 2] > 0).astype(np.int64)
    w = np.array([0., 0., 1.])
    same = np.cross(w, g.xyz)
    opp = np.where(pid[:, None] == 1, same, -same)
    c1 = np.array(tp.coherence(g, pid, same, np.random.default_rng(0), n_ref=60, n_cmp=3000))
    c2 = np.array(tp.coherence(g, pid, opp, np.random.default_rng(0), n_ref=60, n_cmp=3000))
    # one rigid rotation: velocities align at short range and turn apart with distance, so C(d) is ~1 only near
    k = c1 != 0
    assert k.sum() >= 3
    assert c1[1] > .95 and np.all(np.diff(c1[k]) < .02)
    assert np.allclose(c2[k], -c1[k])                                 # the same motion reversed on one plate
