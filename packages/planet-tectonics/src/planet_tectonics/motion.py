"""PlateSim, motion: boundary classification, ocean-ocean polarity per plate pair, the force law and plate motion."""
import numpy as np

from planet_tectonics.boundaries import classify_boundaries, classify_boundaries_ring
from planet_tectonics.params import COLLISION, RIDGE, SUBDUCTION


class MotionMixin:
    """PlateSim, motion: boundary classification, ocean-ocean polarity per plate pair, the force law and plate motion."""

    def _classify(self):
        g = self.g
        K = len(self.plates)
        pid = np.array([pl.id for pl in self.plates])
        slot = np.zeros(pid.max()+1, np.int64); slot[pid] = np.arange(K)
        omega = np.stack([pl.omega(g.radius) for pl in self.plates])
        share = np.array([(g.area*(pl.mask & pl.cont)).sum() for pl in self.plates])/g.area.sum()
        age = self.t-self.world['born']
        score = self._score(self.world['cont'], age, share[slot[self.owner]])
        # owner may refer to plates removed this step only if they own nothing, so slot lookup is safe
        bn = self.p['boundary_normal']
        if bn not in ('node', 'ring'):
            raise ValueError(f"boundary_normal {bn!r}: 'node' or 'ring'")
        fn = classify_boundaries_ring if bn == 'ring' else classify_boundaries
        bt, rate, trench = fn(self.owner, slot, omega, g.xyz, g.nbr, score, self.world['cont'], g.radius,
                              self.p['boundary_min_mm_yr']*1e-3)
        tp = self.p['trench_polarity']
        if tp == 'pair':
            trench = self._pair_polarity(bt, trench)
        elif tp != 'node':
            raise ValueError(f"trench_polarity {tp!r}: 'node' or 'pair'")
        self.btype, self.conv_rate, self.trench = bt, rate, trench


    def _pair_polarity(self, bt, trench):
        """trench_polarity 'pair' (B-d2): oceanic subduction nodes facing an oceanic node of another plate take the
        polarity the overlap resolution used for that pair this step (the overriding plate's nodes are not trench nodes,
        the other plate's are); pairs that did not overlap keep the node rule."""
        g = self.g
        now = getattr(self, 'polarity_now', {})
        if not now:
            return trench
        own, cont = self.owner, self.world['cont']
        I = np.flatnonzero((bt == SUBDUCTION) & ~cont)
        if I.size == 0:
            return trench
        nb = g.nbr[I]
        nbc = np.clip(nb, 0, None)
        oth = (nb >= 0) & (own[nbc] != own[I, None]) & ~cont[nbc]
        has = oth.any(1)
        J = nbc[np.arange(I.size), np.argmax(oth, 1)]
        I, J = I[has], J[has]
        a, b = own[I], own[J]
        win = np.array([now.get((int(x), int(y)), -1) for x, y in zip(np.minimum(a, b), np.maximum(a, b))], np.int64)
        m = win >= 0
        out = trench.copy()
        out[I[m]] = win[m] != a[m]
        self.erosion_totals['polarity_flipped_nodes'] = self.erosion_totals.get('polarity_flipped_nodes', 0)+int(
            (out[I[m]] != trench[I[m]]).sum())
        self.erosion_totals['polarity_pair_nodes'] = self.erosion_totals.get('polarity_pair_nodes', 0)+int(m.sum())
        return out


    def _driving(self, own, speed, axis, relax=None):
        """The usual forces on the plate made of the world nodes own, one step from its speed and axis: the speed
        relaxes toward the target set by slab pull, ridge push, continental drag and collision resistance, and the pole
        turns so the plate moves toward its trenches (slab pull) and away from its ridges (ridge push). Returns the new
        speed and axis."""
        g, p = self.g, self.p
        dt = p['dt_myr']
        h = np.sqrt(g.area.mean())
        c = (g.area[own, None]*g.xyz[own]).sum(0); c /= np.linalg.norm(c)
        tr = own & self.trench
        ri = own & (self.btype == RIDGE)
        coll = own & (self.btype == COLLISION)
        L_slab = tr.sum()*h/2e3                          # boundary nodes sit on both sides of a line: halve
        L_ridge = ri.sum()*h/2e3
        L_coll = coll.sum()*h/2e3
        cf = (g.area*(own & self.world['cont'])).sum()/(g.area*own).sum()
        target = (p['v_min_cm_yr']+(p['v_slab_cm_yr']-p['v_min_cm_yr'])*min(1., L_slab/p['slab_ref_km'])
                  + p['v_ridge_cm_yr']*min(1., L_ridge/p['ridge_ref_km']))
        target *= (1-p['cont_drag']*cf)/(1+L_coll/p['collision_ref_km'])
        speed = speed+(target-speed)*(min(1., dt/p['tau_speed_myr']) if relax is None else relax)
        pull = np.cross(c, g.xyz[tr]).sum(0)+p['ridge_torque']*np.cross(g.xyz[ri], c).sum(0)
        nrm = np.linalg.norm(pull)
        if nrm > 0:
            b = min(1., dt/p['tau_axis_myr']) if relax is None else relax
            a = (1-b)*axis+b*pull/nrm
            axis = a/np.linalg.norm(a)
        return speed, axis


    def _update_motion(self):
        for pl in self.plates:
            own = self.owner == pl.id
            if not own.any():
                continue
            pl.speed, pl.axis = self._driving(own, pl.speed, pl.axis)


    # ------------------------------------------------------------------ plate splitting
    def _omega_vec(self, pl):
        return pl.axis*pl.speed


    def _set_omega(self, pl, v):
        sp = float(np.linalg.norm(v))
        pl.speed = max(sp, 1e-3)
        if sp > 1e-9:
            pl.axis = v/sp
