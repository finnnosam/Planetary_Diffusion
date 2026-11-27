"""Pass 51: E3's record (as precheck-p50-geoid/e3_loads.py) with each parcel's subducted plate area carried
(summed when lower-mantle parcels merge; positions, depths and birth times stay mass-weighted, so the mass record is
identical to D's). Saves the present-day load layers of mass (no alpha) and of area, and mass / area by layer.
argv: trenches.npz ages.npz out.npz [lm_sink_mm_yr]"""
import sys
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import mantle_flow as mf
from planet_grid import fields
from planet_tectonics.params import PARAMS as P

R = 6371e3


class AreaRecord(mf.SlabRecord):
    def __init__(self, *a):
        super().__init__(*a); self.area = np.zeros(0)

    def deposit(self, xyz, dirn, mass, speed, t, area=None):
        keep = mass > 0
        super().deposit(xyz, dirn, mass, speed, t)
        if keep.any():
            self.area = np.concatenate([self.area, area[keep]])

    def _take(self, m):
        super()._take(m); self.area = self.area[m]

    def merge(self, graph, z_cut):
        lm = self.z >= mf.Z_660_M
        if lm.sum() < 2:
            return
        I = np.flatnonzero(lm)
        node = fields.nearest_nodes(graph, self.pos[I])
        lay = np.floor((self.z[I]-z_cut)/mf.LAYER_M).astype(np.int64)
        u, inv = np.unique(node.astype(np.int64)*4096+lay, return_inverse=True)
        if u.size == I.size:
            return
        A = np.bincount(inv, self.area[I], u.size)
        area = np.concatenate([self.area[~lm], A])
        super().merge(graph, z_cut)
        assert area.size == self.mass.size
        self.area = area


T = np.load(sys.argv[1]); Ag = np.load(sys.argv[2]); out = sys.argv[3]
sink = float(sys.argv[4]) if len(sys.argv) > 4 else P['lm_sink_mm_yr']
age = Ag['age_myr']; assert age.shape == T['t'].shape
g = CubeSphere(217)
z_cut = P['flow_z_cut_km']*1e3; lmax = int(P['flow_lmax'])
radii = mf.load_radii(R, z_cut)
rec = AreaRecord(P['slab_dip_deg'], sink)
tt = T['t']
for t in range(250, -1, -1):
    rec.advance(1., R)
    m = (tt == t) & np.isfinite(age) & (T['conv_mm_yr'] > 0)
    conv = T['conv_mm_yr'][m]*1e-3
    area = conv*1e6*T['length_km'][m]*1e3
    mass = mf.thermal_mass_excess(age[m], P['rho_mantle'], P['rho_water'])*area
    rec.deposit(T['xyz'][m], T['n_over'][m], mass, conv, -float(t), area)
    rec.merge(g, z_cut)
    if t % 25 == 0:
        print(t, rec.size, rec.total_mass(), rec.area.sum(), flush=True)
mass_loads = rec.layer_loads(z_cut, len(radii), R, lmax, 1, None)
keep = rec.mass.copy(); rec.mass = rec.area.copy()
area_loads = rec.layer_loads(z_cut, len(radii), R, lmax, 1, None)
area_by_layer = rec.mass_by_layer(z_cut, len(radii)); rec.mass = keep
np.savez_compressed(out, radii=radii, lmax=lmax, sink=sink, loads_0=mass_loads, area_loads_0=area_loads,
                    mass_by_layer_0=rec.mass_by_layer(z_cut, len(radii)), area_by_layer_0=area_by_layer,
                    held=rec.total_mass(), deposited=rec.deposited_kg, exited=rec.exited_kg)
print('done', rec.total_mass(), rec.area.sum())
