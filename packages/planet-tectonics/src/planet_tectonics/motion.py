"""PlateSim, motion: boundary classification, ocean-ocean polarity per plate pair, the force law and plate motion."""
import numpy as np

from planet_grid import fields

from planet_tectonics.boundaries import boundary_frame_ring, classify_boundaries, classify_boundaries_ring
from planet_tectonics.params import COLLISION, CONVERGENT, RIDGE, SUBDUCTION

BOUNDARY_NODE_H = .796             # km of boundary line per boundary node on one side, in units of the node spacing h
                                   # (8-neighbour cube-sphere; 24 random great circles at n = 96 and 217: median 0.796,
                                   # p10-p90 0.72-0.92; tests/test_plates.py measures it)
YR_S = 3.15576e7
SLAB_BAND_RINGS = 2                # N: the slab count is kept this many rings into the overriding plate, so it survives
                                   # the edge moving by up to two cells (~92 km at n = 217) in a step


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
        si = self.p['subduction_initiation']
        if si == 'convergence':
            bt, trench = self._slab_gate(bt, rate, trench)
        elif si != 'instant':
            raise ValueError(f"subduction_initiation {si!r}: 'instant' or 'convergence'")
        self.btype, self.conv_rate, self.trench = bt, rate, trench


    def _spread_same_plate(self, v, rings):
        """Max of v over each node and its neighbours on the same plate, repeated rings times (v >= 0)."""
        g, own = self.g, self.owner
        nb = np.where(g.nbr >= 0, g.nbr, np.arange(g.size)[:, None])
        same = own[nb] == own[:, None]
        for _ in range(rings):
            v = np.maximum(v, np.where(same, v[nb], 0).max(1))
        return v


    def _slab_gate(self, bt, rate, trench):
        """subduction_initiation 'convergence' (B-d3): a subduction boundary has a slab only where its overriding side
        has counted slab_init_km of convergence (Gurnis et al. 2004); elsewhere it is CONVERGENT (no slab, no trench).
        The count lives on the overriding plate's raster ('slab', km) in a band SLAB_BAND_RINGS deep behind its edge:
        each overriding node reads the band's largest value within that reach on its plate (so the count follows the
        edge, and extends along it by up to two reaches per step as the boundary lengthens: lateral propagation of a
        slab), adds its closing rate x dt, and writes the band back. An
        edge typed ridge clears the band behind it (the slab detaches where the boundary opens), except where the band
        also lies behind a convergent edge. A pair made by the subduction-initiation event starts with its slab."""
        g, p = self.g, self.p
        n = g.size
        own = self.owner
        R = SLAB_BAND_RINGS
        D = p['slab_init_km']
        sub = bt == SUBDUCTION
        over = sub & ~trench
        band = self._spread_same_plate(over.astype(np.float64), R) > 0
        ridge_band = self._spread_same_plate((bt == RIDGE).astype(np.float64), R) > 0
        clear = ridge_band & ~band
        look = band | clear
        byid = {pl.id: pl for pl in self.plates}
        ids = np.flatnonzero(look)
        rows = {}
        V = np.zeros(n)
        for k in np.unique(own[ids]):
            pl = byid.get(int(k))
            if pl is None:
                continue
            I = ids[own[ids] == k]
            r = fields.nearest_nodes(g, g.xyz[I]@pl.Q)
            rows[int(k)] = (pl, I, r)
            V[I] = pl.slab[r]
            c = clear[I]
            pl.slab[r[c]] = 0
        V[clear] = 0
        S = self._spread_same_plate(np.where(band, V, 0.), R)
        S = np.where(over, np.minimum(S+np.maximum(rate, 0.)*1e3*p['dt_myr'], 10*D), 0.)   # m/yr x 1e3 x Myr = km
        seeds = getattr(self, 'slab_seed_pairs', set())
        if seeds:                                        # spontaneous initiation: the new pair starts with its slab
            I = np.flatnonzero(over)
            nb = g.nbr[I]
            nbc = np.clip(nb, 0, None)
            oth = (nb >= 0) & (own[nbc] != own[I, None])
            J = nbc[np.arange(I.size), np.argmax(oth, 1)]
            hit = oth.any(1) & np.array([(min(a, b), max(a, b)) in seeds for a, b in zip(own[I], own[J])], bool)
            S[I[hit]] = np.maximum(S[I[hit]], D)
            done = {(min(a, b), max(a, b)) for a, b in zip(own[I[hit]], own[J[hit]])}
            self.slab_seed_pairs = seeds-done
        W = self._spread_same_plate(S, R)
        for pl, I, r in rows.values():
            m = band[I] & pl.mask[r]                     # only onto the plate's own crust
            np.maximum.at(pl.slab, r[m], W[I[m]].astype(np.float32))
        has = over & (S >= D)
        near = has.copy()
        nb = np.where(g.nbr >= 0, g.nbr, np.arange(n)[:, None])
        for _ in range(R):                              # the downgoing side within reach of an overriding node with a slab
            near = near | near[nb].any(1)
        slab = has | (sub & trench & near)
        out_bt = np.where(sub & ~slab, CONVERGENT, bt).astype(bt.dtype)
        out_tr = trench & slab
        ET = self.erosion_totals
        h = np.sqrt(g.area.mean())/1e3
        ET['convergent_no_slab_node_km_steps'] = ET.get('convergent_no_slab_node_km_steps', 0.)+float((out_bt == CONVERGENT).sum()*h/2)
        return out_bt, out_tr


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


    def _torque_omega(self, plates):
        """force_law 'torque' (forty-third pass): the angular velocities (rad/yr) that balance the torques on all
        plates at once (Forsyth & Uyeda 1975; no inertia), SI inside:
          sum over the plate's trench nodes of F_sp L R (x x n)  [slab pull, n the normal out of the plate]
        + sum over its overriding nodes of F_su L R (x x n)      [slab suction toward the trench]
        + sum over its oceanic ridge nodes of F_rp L R (x x -n)  [ridge push]
        = C_b R^2 sum_nodes a (1 + (keel-1) cont) (I - x x^T) w   [basal drag; keels under continental nodes]
        + C_s L R^2 sum_trench q q^T w                           [slab resistance, q = x x n]
        + C_c L R^2 sum_{collision, convergent} q q^T (w - w_other)   [convergence without a slab]
        L = BOUNDARY_NODE_H h per boundary node. Coefficients from the fit to Earth's plates on Bird's boundaries
        (planet_reviews.force_balance_earth)."""
        g, p = self.g, self.p
        K = len(plates)
        pid = np.array([pl.id for pl in plates])
        slot = np.full(max(pid.max(), int(self.owner.max()))+1, -1, np.int64); slot[pid] = np.arange(K)
        omega = np.stack([pl.omega(g.radius) for pl in plates])
        full = np.zeros(slot.size, np.int64); full[pid] = np.arange(K)
        other, nrm = boundary_frame_ring(self.owner, full, omega, g.xyz, g.nbr, g.radius)
        R = g.radius
        Ln = BOUNDARY_NODE_H*np.sqrt(g.area.mean())
        x = g.xyz
        row = slot[self.owner]
        cont = self.world['cont']
        w = g.area*(1+(p['keel_drag_factor']-1)*cont)
        A = np.zeros((K, K, 3, 3)); rhs = np.zeros((K, 3))
        ok = row >= 0
        sw = np.bincount(row[ok], w[ok], K)
        for i in range(3):
            for j in range(3):
                A[np.arange(K), np.arange(K), i, j] -= p['basal_drag_pa_s_m']*R**2*np.bincount(row[ok], (w*x[:, i]*x[:, j])[ok], K)
            A[np.arange(K), np.arange(K), i, i] += p['basal_drag_pa_s_m']*R**2*sw
        Fsp = p['slab_pull_n_m']
        Cs = Fsp/(p['slab_terminal_cm_yr']*1e-2/YR_S)
        Fsu = p['slab_suction_frac']*Fsp
        Cc = p['convergent_coupling']*Cs
        q = np.cross(x, nrm)
        bt = self.btype
        def add_force(m, F, sign):
            if F:
                np.add.at(rhs, row[m], sign*F*Ln*R*q[m])
        def add_drag(m, C, couple):
            if not C or not m.any():
                return
            qq = C*Ln*R**2*q[m, :, None]*q[m, None, :]
            np.add.at(A, (row[m], row[m]), qq)
            if couple:
                np.add.at(A, (row[m], slot[other[m]]), -qq)
        b = (other >= 0) & ok
        tr = b & self.trench
        add_force(tr, Fsp, 1.); add_drag(tr, Cs, False)
        add_force(b & (bt == SUBDUCTION) & ~self.trench, Fsu, 1.)
        add_force(b & (bt == RIDGE) & ~cont, p['ridge_push_n_m'], -1.)
        cv = b & ((bt == COLLISION) | (bt == CONVERGENT))
        cv &= slot[np.where(other >= 0, other, 0)] >= 0
        add_drag(cv, Cc, True)
        M = A.transpose(0, 2, 1, 3).reshape(3*K, 3*K)
        try:
            w = np.linalg.solve(M, rhs.reshape(-1))
        except np.linalg.LinAlgError:                    # a plate of one or two nodes has no drag against spin about
            w = np.linalg.lstsq(M, rhs.reshape(-1), rcond=None)[0]     # its own axis: the least-norm balance
        return w.reshape(K, 3)*YR_S                      # rad/s -> rad/yr


    def _update_motion(self):
        fl = self.p['force_law']
        if fl == 'torque':
            plates = [pl for pl in self.plates if (self.owner == pl.id).any()]
            if not plates:
                return
            wt = self._torque_omega(plates)
            b = min(1., self.p['dt_myr']/self.p['tau_speed_myr'])
            for pl, w in zip(plates, wt):
                self._set_omega(pl, (1-b)*self._omega_vec(pl)+b*w*self.g.radius*1e2)    # speed x axis in cm/yr
            return
        if fl != 'lengths':
            raise ValueError(f"force_law {fl!r}: 'lengths' or 'torque'")
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
