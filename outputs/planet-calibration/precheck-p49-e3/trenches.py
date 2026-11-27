"""E3 step 1 (pass 49): Earth's subduction zones through time from Müller et al. (2019) v3.0 (v2.02), via pygplates.

Per Myr t (0 ... 250 Ma): resolve the topologies (boundaries and deforming networks), walk every SubductionZone shared
sub-segment, take its overriding and subducting plates (pygplates' own polarity), tessellate it at STEP_DEG, and at
each arc midpoint record the position (absolute frame = the model's mantle frame), the arc length, the unit normal
toward the overriding side, and the convergence rate normal to the trench, (v_sub - v_over) . n_over, from 1 Myr
stage rotations (t+1 -> t) of the two plates' reconstruction plate ids.
Writes trenches.npz: t, xyz, length_km, n_over (xyz), conv_mm_yr, sub_pid, over_pid.
"""
import sys
import numpy as np
import pygplates

STEP_DEG = 0.5
R_KM = 6371.0


def run(rot_file, topo_file, times, out):
    rm = pygplates.RotationModel(rot_file)
    feats = pygplates.FeatureCollection(topo_file)
    cols = {k: [] for k in ('t', 'xyz', 'length_km', 'n_over', 'conv_mm_yr', 'sub_pid', 'over_pid')}
    stats = []
    for t in times:
        resolved, shared = [], []
        pygplates.resolve_topologies(feats, rm, resolved, float(t), shared)
        L_all = 0.; L_used = 0.; nskip = 0
        for sec in shared:
            if sec.get_feature().get_feature_type() != pygplates.FeatureType.gpml_subduction_zone:
                continue
            for ss in sec.get_shared_sub_segments():
                geom = ss.get_resolved_geometry()
                L = geom.get_arc_length()*R_KM
                L_all += L
                res = ss.get_overriding_and_subducting_plates(True)
                if res is None:
                    nskip += 1
                    continue
                over, sub, pol = res
                po = over.get_feature().get_reconstruction_plate_id()
                ps = sub.get_feature().get_reconstruction_plate_id()
                tess = geom.to_tessellated(np.radians(STEP_DEG))
                pts = np.array([p.to_xyz() for p in tess.get_points()])
                if len(pts) < 2:
                    continue
                a, b = pts[:-1], pts[1:]
                mid = a+b; mid /= np.linalg.norm(mid, axis=1, keepdims=True)
                gc = np.cross(a, b); ln = np.linalg.norm(gc, axis=1)
                ok = ln > 0
                a, b, mid, gc, ln = a[ok], b[ok], mid[ok], gc[ok], ln[ok]
                left = gc/ln[:, None]                          # points to the left of travel a -> b
                n_over = left if pol == 'Left' else -left
                arc_km = np.arcsin(np.clip(ln, 0, 1))*R_KM
                vel = {}
                for pid in (po, ps):
                    stage = rm.get_rotation(float(t), pid, float(t)+1.)
                    v = pygplates.calculate_velocities([pygplates.PointOnSphere(*m) for m in mid], stage, 1.,
                                                       pygplates.VelocityUnits.kms_per_my)
                    vel[pid] = np.array([x.to_xyz() for x in v])   # km/Myr = mm/yr
                conv = np.einsum('ij,ij->i', vel[ps]-vel[po], n_over)
                n = len(mid)
                cols['t'].append(np.full(n, float(t))); cols['xyz'].append(mid); cols['length_km'].append(arc_km)
                cols['n_over'].append(n_over); cols['conv_mm_yr'].append(conv)
                cols['sub_pid'].append(np.full(n, ps)); cols['over_pid'].append(np.full(n, po))
                L_used += arc_km.sum()
        c = np.concatenate(cols['conv_mm_yr'][-0:]) if False else None
        stats.append((t, L_all, L_used, nskip))
        print(f't={t:5.0f} Ma: subduction length {L_all/1e3:.1f}k km, with plates {L_used/1e3:.1f}k km, sub-segments '
              f'without plates {nskip}', flush=True)
    np.savez_compressed(out, **{k: np.concatenate(v) for k, v in cols.items()}, stats=np.array(stats))


if __name__ == '__main__':
    rot, topo, out = sys.argv[1:4]
    t0, t1 = (int(x) for x in sys.argv[4:6]) if len(sys.argv) > 5 else (0, 250)
    run(rot, topo, range(t0, t1+1), out)
