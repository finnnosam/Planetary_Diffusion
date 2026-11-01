"""D, the slab-driven mantle flow (planet_tectonics.mantle_flow; journal 66 section 5, tests U1-U6)."""
import numpy as np
import pytest

from planet_grid.paths import DATA
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_tectonics import plates as P

R = 6371e3


def test_u1_eigen_propagator_matches_rk4():
    for l in (1, 2, 8, 31):
        for e in (.01, 1., 100.):
            x = np.log(6371./5700.)
            A, B = mf.propagator_eig(l, e, x), mf.propagator_rk4(l, e, x)
            assert np.abs(A-B).max() <= 1e-9*max(1., np.abs(B).max())
    w = np.sort(np.linalg.eigvals(mf.system_matrix(5, 3.)).real)
    assert np.allclose(w, [-7, -5, 4, 6])                      # -l-2, -l, l-1, l+1


def _iso(lmax, depths_km):
    return mf.FlowKernels(R, [R-z*1e3 for z in depths_km], np.array([mf.CMB_RADIUS_M/R]), np.array([1e21]), lmax)


def test_u2_kernel_limits_sign_and_decay():
    K = _iso(30, [1., 100., 1500., 2500., 2880.])
    kh = K.KH
    assert (kh[1:, :4] < 0).all()                              # a dense sheet pulls the surface down
    assert abs(kh[2, 0]+1) < .02                               # a load at the surface: isostatic, |K| -> 1
    assert abs(kh[2, 1]+.96) < .02 and abs(kh[2, 2]+.40) < .02 and abs(kh[2, 3]+.09) < .02   # isoviscous l = 2
    assert abs(kh[2, 4]) < .02                                 # at the CMB: compensated there, -> 0
    assert abs(kh[30, 2]) < abs(kh[4, 2]) < abs(kh[2, 2])      # faster decay with depth at high degree
    assert abs(kh[30, 2]) < 3e-3


def test_u3_gradient_matches_finite_differences():
    lmax = 12
    rng = np.random.default_rng(1)
    a = (rng.normal(size=(lmax+1)*(lmax+2)//2)+1j*rng.normal(size=(lmax+1)*(lmax+2)//2))
    a[mf.alm_m(lmax) == 0] = a[mf.alm_m(lmax) == 0].real
    x = rng.normal(size=(50, 3)); x /= np.linalg.norm(x, axis=1)[:, None]
    gr = mf.gradient(a, lmax, x)
    assert np.abs((gr*x).sum(1)).max() < 1e-9                  # tangent
    h = 1e-5
    for k in range(3):
        e = np.zeros(3); e[k] = 1.
        t = e-(x@e)[:, None]*x                                 # tangent direction
        tn = t/np.linalg.norm(t, axis=1)[:, None]
        xp = x+h*tn; xp /= np.linalg.norm(xp, axis=1)[:, None]
        xm = x-h*tn; xm /= np.linalg.norm(xm, axis=1)[:, None]
        fd = (mf.synthesis(a, lmax, xp)-mf.synthesis(a, lmax, xm))/(2*h)
        assert np.abs(fd-(gr*tn).sum(1)).max() < 1e-6*np.abs(fd).max()
    # the gradient power equals the mean square of the synthesised gradient (quadrature on many random points)
    xx = rng.normal(size=(400000, 3)); xx /= np.linalg.norm(xx, axis=1)[:, None]
    ms = (mf.gradient(a, lmax, xx)**2).sum(1).mean()
    assert abs(ms/mf.gradient_power(a, lmax)-1) < .02


def test_u4_record_conserves_mass_and_sinks_at_the_set_rates():
    g = CubeSphere(24)
    rec = mf.SlabRecord(58., 12.)
    x = np.array([[1., 0, 0], [0, 1., 0], [1., 0, 0]]); d = np.array([[0, 1., 0], [0, 0, 1.], [0, 1., 0]])
    rec.deposit(x, d, np.array([1e15, 2e15, 0.]), np.array([.05, .02, .05]), -100.)
    assert rec.size == 2 and rec.deposited_kg == 3e15            # a parcel without mass is not kept
    rec.advance(1., R)
    sd, cd = np.sin(np.radians(58)), np.cos(np.radians(58))
    assert np.allclose(rec.z, [5e4*sd, 2e4*sd])
    ang = np.arccos(np.clip((rec.pos*x[:2]).sum(1), -1, 1))*R
    assert np.allclose(ang, [5e4*cd, 2e4*cd], rtol=1e-9)
    assert np.allclose((rec.pos*d[:2]).sum(1) > 0, True)        # toward the overriding side
    # the first parcel reaches 660 km after 660 / (50 sin 58) Myr; the rest of that step it sinks at 12 mm/yr
    t_cross = mf.Z_660_M/(.05*sd)/1e6
    rec.advance(t_cross-1+.5, R)
    assert abs(rec.z[0]-(mf.Z_660_M+.5e6*.012)) < 1e-6
    rec.advance(10., R)
    assert abs(rec.z[0]-(mf.Z_660_M+10.5e6*.012)) < 1e-6
    rec.deposit(rec.pos[:1].copy(), rec.dirn[:1].copy(), np.array([5e14]), np.array([.05]), -80.)
    rec.z[-1] = rec.z[0]+1e3
    rec.merge(g, 300e3)
    assert abs(rec.total_mass()+rec.exited_kg-rec.deposited_kg) < 1e-6*rec.deposited_kg
    lm = rec.z >= mf.Z_660_M
    assert lm.sum() == 1 and abs(rec.mass[lm][0]-1.5e15) < 1 and abs(rec.born[lm][0]-(-100*2+-80)/3*1.) < 1e-9
    rec.advance(300., R)                                         # everything leaves at 2,850 km
    assert rec.size == 0 and abs(rec.exited_kg-rec.deposited_kg) < 1e-6*rec.deposited_kg


def test_u4_layer_loads_sum_to_the_mass():
    rec = mf.SlabRecord()
    rng = np.random.default_rng(0)
    x = rng.normal(size=(30, 3)); x /= np.linalg.norm(x, axis=1)[:, None]
    rec.deposit(x, np.zeros((30, 3)), np.full(30, 1e16), np.full(30, .05), 0.)
    rec.z[:] = np.linspace(200e3, 2800e3, 30)
    radii = mf.load_radii(R, 300e3)
    L = rec.layer_loads(300e3, len(radii), R, 8)
    # degree 0: Y00 = 1/sqrt(4 pi); sum over layers of a00 sqrt(4 pi) r^2 is the mass below the cut
    m = (L[:, 0].real*np.sqrt(4*np.pi)*radii**2).sum()
    assert abs(m-1e16*(rec.z >= 300e3).sum())/m < 1e-9
    assert abs(rec.mass_by_layer(300e3, len(radii)).sum()-m)/m < 1e-9


def test_u5_random_is_the_default_and_slabs_keeps_the_hotspots():
    assert P.PARAMS['dyn_rule'] == 'random' and P.PARAMS['basal_traction'] is False and P.PARAMS['slab_alpha'] == 'thermo'
    if not (DATA/'mantle'/'visc.sh08').exists():
        pytest.skip('data/mantle/visc.sh08 missing')
    g = CubeSphere(24)
    a = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'erosion': False, 'threads': 1})
    b = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'erosion': False, 'threads': 1, 'dyn_rule': 'slabs'})
    assert a.mantle.flow is None and b.mantle.flow is not None
    assert [h['pos'].tolist() for h in a.mantle.hotspots] == [h['pos'].tolist() for h in b.mantle.hotspots]
    for _ in range(3):
        a.step(); b.step()
    rec = b.mantle.flow.record
    assert rec.deposited_kg > 0 and abs(rec.total_mass()+rec.exited_kg-rec.deposited_kg) < 1e-9*rec.deposited_kg
    st = b.state()
    for k in ('flow_topography_alm_m', 'flow_velocity_alm_m_yr', 'slab_mass_by_layer_kg', 'dynamic_topography_air_m'):
        assert k in st
    assert 'flow_topography_alm_m' not in a.state()


def test_basal_traction_with_a_rigid_flow_moves_the_plates_with_it():
    """No new coefficient: if the flow is one rigid rotation and nothing else drives, every plate takes that rotation."""
    if not (DATA/'mantle'/'visc.sh08').exists():
        pytest.skip('data/mantle/visc.sh08 missing')
    g = CubeSphere(24)
    sim = P.PlateSim(g, seed=3, params={'water_volume_km3': None, 'erosion': False, 'threads': 1, 'dyn_rule': 'slabs',
                                        'force_law': 'torque', 'basal_traction': True, 'slab_pull_n_m': 0.,
                                        'ridge_push_n_m': 0.})
    sim.step()
    w0 = np.array([.2, -.5, 1.])*1e-9                            # rad/yr
    sim.mantle.flow.v_xyz = np.cross(w0, g.xyz)*g.radius
    plates = [pl for pl in sim.plates if (sim.owner == pl.id).any()]
    w = sim._torque_omega(plates)
    assert np.allclose(w, w0[None, :], rtol=1e-6, atol=1e-15)


@pytest.mark.skipif(not (DATA/'mantle'/'smean2.31.m.ab').exists(), reason='data/mantle missing')
def test_u6_earth_precheck_numbers_reproduce():
    from planet_reviews import mantle_flow_earth as mfe
    s = mfe.summary(*mfe.earth_flow('smean2', True, 300.))
    assert abs(s['rms_m']/222-1) < .01 and abs(s['flow_rms_cm_yr']/1.35-1) < .01
    assert abs(s['share_l1_3']-.84) < .01


def test_thermo_ratio_is_the_pass51_table():
    z = mf.THERMO_Z0_KM+mf.THERMO_DZ_KM*np.arange(len(mf.THERMO_F))
    assert len(mf.THERMO_F) == len(mf.load_radii(R, 300e3)) == 51 and z[-1] == 2825.
    assert np.allclose((R-mf.load_radii(R, 300e3))/1e3, z)          # tabulated at the load-layer mid-depths
    f = mf.thermo_ratio(z)
    assert np.allclose(f, mf.THERMO_F) and abs(mf.thermo_ratio(325.)-.7106) < 1e-12
    lm = z > 700.
    assert (np.diff(f[lm]) < 0).all()                               # falls through the lower mantle
    r = f/mf.alpha_ratio(z)                                         # above the ambient-alpha factor below 800 km,
    assert (r[z > 800.] > 1.05).all() and (np.diff(r[lm]) > 0).all() and 1.55 < r[-1] < 1.65   # growing with depth
    assert mf.thermo_ratio(100.) == f[0] and mf.thermo_ratio(2891.) == f[-1]


def test_alpha_ratio_follows_its_sources():
    a = mf.alpha_ratio(np.array([0., 400., 660., 670., 2891.]))*mf.ALPHA_GDH1
    assert np.allclose(a, [3.1e-5, 2.378e-5, 2.079e-5, 2.5e-5, 1.0e-5], rtol=1e-9)
    z = np.linspace(670., 2891., 50)
    assert (np.diff(mf.alpha_ratio(z)) < 0).all()                   # falls monotonically through the lower mantle
    rec = mf.SlabRecord()
    rec.deposit(np.array([[0., 0, 1.]]), np.zeros((1, 3)), np.array([1e16]), np.array([.05]), 0.)
    rec.z[:] = 1520e3
    radii = mf.load_radii(R, 300e3)
    f = mf.alpha_ratio((R-radii)/1e3)
    L0, L1 = rec.layer_loads(300e3, len(radii), R, 4), rec.layer_loads(300e3, len(radii), R, 4, factor=f)
    j = int((1520e3-300e3)//mf.LAYER_M)
    assert np.allclose(L1[j], L0[j]*f[j]) and abs(f[j]-mf.alpha_ratio(1525.)) < 1e-12
