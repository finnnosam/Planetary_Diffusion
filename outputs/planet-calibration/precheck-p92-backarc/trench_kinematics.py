"""Pass 92 (b): the upper plate's regime on the planet's continental overriding trench under three trench laws, read
from saved snapshots (plate rotations in t*-plates.json, boundary types, crust): the overriding side is a continental
subduction-boundary node with an oceanic neighbour on another plate; u its plate's speed toward that neighbour (the
trench normal as PlateSim._trenchward_speed), Vsub the subducting plate's speed toward the upper plate (landward).
Laws: rollback 17 mm/yr (Oncken 2006), 9 mm/yr (Lallemand 2008 Table 2 mean, SB04), and Lallemand 2008's neutral line
in SB04, rollback = 14 - 0.24 Vsub (mm/yr). Shares of trench nodes shortening / extending beyond 5 mm/yr against
Earth's continental upper plates (Lallemand 2008 Table 1: compression 0.21, extension 0.11). Args: arm folders.
Run from the repo root."""
import json, sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
g = CubeSphere(217)
R = g.radius
SUB = 2
LAWS = (('rollback 17', lambda vs: 17.), ('rollback 9', lambda vs: 9.), ('Lallemand SB04', lambda vs: 14.-.24*vs))
for arm in sys.argv[1:]:
    acc = {k: [] for k, _ in LAWS}; U, VS = [], []
    for sd in sorted(Path(arm).glob('seed*')):
        for f in sorted((sd/'snapshots').glob('t*-plates.json')):
            q = np.load(str(f).replace('-plates.json', '.npz'))
            if float(q['time_myr']) < -40.5:
                continue
            om = {r['id']: np.asarray(r['axis'])*r['speed_cm_yr']*1e-2/R for r in json.loads(f.read_text())}
            # speed_cm_yr is the plate's speed; the saved axis is in the world frame at the save
            own = q['plate_id'].astype(np.int64); cont = q['continental'].astype(bool); bt = q['boundary_type']
            idx = np.flatnonzero((bt == SUB) & cont)
            nb = g.nbr[idx]; nbs = np.maximum(nb, 0)
            other = (nb >= 0) & (own[nbs] != own[idx][:, None])
            ocean = other & ~cont[nbs]
            keep = ocean.any(1)
            idx, nb, nbs, other, ocean = idx[keep], nb[keep], nbs[keep], other[keep], ocean[keep]
            x = g.xyz[idx]
            disp = np.where(other[..., None], g.xyz[nbs]-x[:, None], 0.).sum(1)
            disp -= np.einsum('ij,ij->i', disp, x)[:, None]*x
            n = disp/np.maximum(np.linalg.norm(disp, axis=1, keepdims=True), 1e-15)
            w = np.stack([om.get(int(k), np.zeros(3)) for k in own[idx]])
            u = np.einsum('ij,ij->i', np.cross(w, x)*R, n)*1e3                       # mm/yr toward the trench
            vs = []
            for r_, (row, oc) in enumerate(zip(nbs, ocean)):
                ws = [om.get(int(own[j]), np.zeros(3)) for j in row[oc]]
                vs.append(np.mean([-np.dot(np.cross(wj, x[r_])*R, n[r_]) for wj in ws])*1e3)   # mm/yr landward
            vs = np.array(vs)
            U.append(u); VS.append(vs)
            for k, law in LAWS:
                acc[k].append(u-law(vs))
    U = np.concatenate(U); VS = np.concatenate(VS)
    print(f'\n{Path(arm).name}: {U.size} continental overriding trench node-snapshots (-40 ... 0 Myr)')
    print(f'  u (upper plate toward the trench) p10 / p50 / p90: {np.percentile(U, 10):.0f} / {np.median(U):.0f} / '
          f'{np.percentile(U, 90):.0f} mm/yr; Vsub (landward) {np.percentile(VS, 10):.0f} / {np.median(VS):.0f} / {np.percentile(VS, 90):.0f}')
    for k, _ in LAWS:
        s = np.concatenate(acc[k])
        print(f'  {k:15s}: shortening > 5 mm/yr {np.mean(s > 5):.2f}, neutral {np.mean(np.abs(s) <= 5):.2f}, '
              f'extending > 5 mm/yr {np.mean(s < -5):.2f}; extension rate p50 / p90 where > 5: '
              f'{np.median(-s[s < -5]) if (s < -5).any() else 0:.0f} / {np.percentile(-s[s < -5], 90) if (s < -5).any() else 0:.0f} mm/yr')
print('\nEarth, continental upper plates (Lallemand 2008 Table 1, 126 transects): compression 0.21, neutral 0.67, extension 0.11')
