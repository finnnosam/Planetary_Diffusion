"""Pass 65 (journal 87): K, the Earth breakup rule (`split_rule` 'earth', `breakup.EARTH_RIFT_PRODUCTION`), and the
margin shortening with the belt at the trench flux (`shortening_match`). Both default since 6 Nov (journal 88)."""
import numpy as np

from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import plates as P
from planet_tectonics.breakup import EARTH_RIFT_PRODUCTION, production_per_mkm2
from planet_tectonics.crust_kernels import margin_shortening


def test_new_rules_are_the_defaults_since_6_nov():
    assert P.PARAMS['split_rule'] == 'earth'           # adopted 6 Nov (journal 88)
    assert P.PARAMS['shortening_match'] is True


def test_earth_production_table_and_its_time_mapping():
    assert len(EARTH_RIFT_PRODUCTION) == 25 and EARTH_RIFT_PRODUCTION[0][0] == 0 and EARTH_RIFT_PRODUCTION[-1][0] == 240
    per = [r[1]/r[2] for r in EARTH_RIFT_PRODUCTION]
    assert abs(np.mean([r[1] for r in EARTH_RIFT_PRODUCTION])-1553) < 1      # the measured mean, km/Myr
    assert production_per_mkm2(0.5) == per[24]          # just after the start: Earth's (240, 250] Ma
    assert production_per_mkm2(90.) == per[15]          # 160 Ma: (150, 160] Ma, the breakup peak
    assert production_per_mkm2(249.5) == per[0] and production_per_mkm2(400.) == per[0] and production_per_mkm2(-3.) == per[24]
    assert max(per) == per[15]


def _segments():
    """Two trench segments: segment 0's belt is missing (one continental node at 300 km, in the far tail of the belt
    Gaussian; an oceanic forearc in front of it), segment 1 has a full belt at 150 km; both have the same foreland."""
    d = np.array([350e3, 350e3, 300e3, 150e3, 150e3, 150e3, 350e3, 350e3])
    near = np.array([0, 0, 0, 1, 1, 1, 1, 1])
    dom = np.ones(8, bool)
    s = np.full(8, .01)
    thick = np.full(8, 40e3)
    area = np.full(8, 2e9)
    return near, d, dom, s, thick, area


def test_match_conserves_volume_and_removes_the_missing_belt_spike():
    near, d, dom, s, thick, area = _segments()
    args = (near, d, dom, s, thick, area, 1e6, 150e3, 60e3, 350e3, 100e3)
    old = margin_shortening(*args)
    new = margin_shortening(*args, match=True)
    for x in (old, new):
        assert abs((x*area).sum()) < 1e-6*np.abs(x*area).sum()
    assert old[2] > 2*old[3]                            # the tail node took most of its foreland's volume
    gb_max = 1/(60e3*np.sqrt(2*np.pi))
    assert new.max() <= .01*1e6*40e3*gb_max*(1+1e-9)    # no node above the trench flux at the belt's peak
    assert new[2] < .2*old[2]                           # segment 0 moves little: its belt is missing
    assert np.allclose(new[3:6], old[3:6], rtol=.5)     # a full belt moves about what it did


def test_match_leaves_a_full_belt_at_the_flux():
    """A segment with a full belt (many nodes at d_belt) and a foreland that can give enough: the belt gets the flux."""
    d = np.r_[np.full(20, 150e3), np.full(20, 350e3)]
    near = np.zeros(40, int); dom = np.ones(40, bool); s = np.full(40, .01); thick = np.full(40, 40e3); area = np.full(40, 2e9)
    new = margin_shortening(near, d, dom, s, thick, area, 1e6, 150e3, 60e3, 350e3, 100e3, match=True)
    gb_max = 1/(60e3*np.sqrt(2*np.pi)); gf_max = 1/(100e3*np.sqrt(2*np.pi))
    assert np.allclose(new[:20], min(1., gf_max/gb_max)*.01*1e6*40e3*gb_max, rtol=.01)   # the smaller side sets it
    assert abs((new*area).sum()) < 1e-6*np.abs(new*area).sum()


def test_earth_split_rule_and_match_run_and_keep_the_ledger():
    g = CubeSphere(48)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'split_rule': 'earth', 'shortening_match': True,
                                        'split_min_area_mkm2': 1.})
    v0 = sim.continental_volume()
    for _ in range(30):
        sim.step()
    assert abs((sim.continental_volume()-v0)-sum(sim.ledger.values())) < 2e-3*v0
    assert sim.t_start == -250.


def test_earth_split_rule_makes_earths_margin_on_average():
    """Expected margin per Myr = production x continental area while the rate per step is small: a quiet bin (Earth's
    (230, 240] Ma, 1.98 km/Myr/Mkm^2), candidate splits drawn but not applied."""
    g = CubeSphere(48)
    # on the as-built start (many continental plates): on the whole start one candidate line through the supercontinent
    # carries ~46,000 km of margin at n = 48, so 300 draws expect only ~2.5 splits (sixty-seventh pass)
    sim = P.PlateSim(g, seed=5, params={'water_volume_km3': None, 'split_rule': 'earth', 'start_plates': 'voronoi',
                                        'start_poles': 'random', 'start_breakup': 'dispersal'})
    sim.step()
    sim.t_start = sim.t-10.5
    prod = production_per_mkm2(10.5)
    last, made = [], []
    real = sim._split_margin_km
    sim._split_margin_km = lambda own, lab: last.append(real(own, lab)) or last[-1]
    # only the rule's own splits: on the coast start oceanic plates hold > 170 Myr floor, and stubbed subduction
    # initiations would fire every pass and be booked at the last candidate's margin (sixty-seventh pass)
    sim.split = lambda pl, part, speed, kind: (made.append(last[-1]) if kind == 'rift' else None) or None
    A = g.area; w = sim.world
    want = sum(prod*A[(sim.owner == pl.id) & w['cont']].sum()/1e12 for pl in sim.plates
               if A[sim.owner == pl.id].sum()/1e12 >= sim.p['split_min_area_mkm2'])
    n = 300
    for _ in range(n):
        sim._split_events_earth()
    got = sum(made)/n
    assert want > 0 and .6*want < got < 1.25*want
