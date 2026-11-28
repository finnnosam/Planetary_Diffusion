"""The geography review (planet_reviews.geography, journal 149): the target measures on a small cube with a known map."""
import json

import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_reviews import geography as geo


def _map(g):
    """A round continent of 50 degrees radius at 600 m, its inner 25 degrees at 100 m, one lake node inside, one
    one-node island in the ocean."""
    c = np.array([0., 0., 1.])
    ang = np.degrees(np.arccos(np.clip(g.xyz@c, -1, 1)))
    z = np.where(ang < 50., 600., -3000.)
    z[ang < 25.] = 100.
    lake = int(np.argmin(np.abs(ang-12.)))
    z[lake] = -50.
    far = np.where(ang > 120.)[0]
    island = int(far[0])
    z[island] = 50.
    return z, lake, island


def test_measures_on_a_known_map():
    g = CubeSphere(48)
    z, lake, island = _map(g)
    m = geo.Geo(g).measures(z)
    land = z > 0
    assert abs(m['land']-g.area[land].sum()/g.area.sum()) < 1e-12
    assert m['holes'] == 1 and m['islands'] == 1
    assert 0 < m['lowflat4'] < m['dry4'] <= 100
    assert m['low200'] > 0 and m['far300'] > 0
    assert 0 < m['low_far300'] <= 100
    # The lake is not the main water body: it makes a coast for far300 but not for far300_sea.
    assert m['far300_sea'] > m['far300']
    water = z <= 0
    assert abs(m['shelf200']-100*g.area[lake]/g.area[water].sum()) < 1e-9
    assert abs(m['coast_step']-3600.) < 1e-9          # most shore edges: 600 m land beside -3,000 m sea
    z2 = z.copy()
    z2[lake] = 100.                                     # fill the lake: no hole, fewer crossings
    m2 = geo.Geo(g).measures(z2)
    assert m2['holes'] == 0 and m2['crossings'] < m['crossings']
    assert abs(m2['far300_sea']-m2['far300']) < 1e-12 and m2['shelf200'] == 0


def test_sweep_folder_reads_its_arms_control_first(tmp_path):
    g = CubeSphere(217)
    z, _, _ = _map(g)
    for arm in ('b-arm', 'a-control'):
        for s in (1, 2):
            d = tmp_path/'sw'/arm/f'seed{s}'/'snapshots'
            d.mkdir(parents=True)
            for t in (-10., -5., 0.):
                np.savez(d/f't{t:+07.1f}.npz', elevation_m=(z+(5. if arm == 'b-arm' else 0.)).astype(np.float32))
    assert [p.split('/')[-1].split('\\')[-1] for p in geo.expand_arms([str(tmp_path/'sw')], 'a-control')] == \
        ['a-control', 'b-arm']
    out = tmp_path/'geo.json'
    geo.main([str(tmp_path/'sw'), '--control', 'a-control', '--window', '10', '--json', str(out)])
    r = json.loads(out.read_text())
    assert set(r) >= {'a-control', 'b-arm'}
    assert set(r['a-control']) == {'seed1', 'seed2'} and r['a-control']['seed1']['n_window'] == 3
