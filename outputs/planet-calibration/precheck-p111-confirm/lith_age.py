"""Pass 111 probe (rule-lith.md): the continental lithosphere age under two arms, every saved snapshot, paired by seed.
Usage (repo root): python lith_age.py <arm A folder> <arm C folder> [seeds, default 1001-1006]
Prints, per snapshot: area-weighted mean lithosphere age of continental nodes (Myr) on A and C, the paired gap C - A
(mean, sd, seeds with C > A), and the area share of continental nodes younger than 200 Myr on each."""
import sys
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere

A, C = Path(sys.argv[1]), Path(sys.argv[2])
seeds = [int(s) for s in sys.argv[3].split(',')] if len(sys.argv) > 3 else list(range(1001, 1007))
area = CubeSphere(217).area
tags = sorted({f.stem for f in (A/f'seed{seeds[0]}'/'snapshots').glob('t*.npz')}, key=lambda t: float(t[1:]))


def read(arm, s, tag):
    d = np.load(arm/f'seed{s}'/'snapshots'/f'{tag}.npz')
    c = d['continental'].astype(bool) & (d['crust_thickness_m'] >= 20e3)
    age = np.nan_to_num(d['lithosphere_age_myr'].astype(float))
    w = area*c
    return float((age*w).sum()/w.sum()), float((w*(age < 200.)).sum()/w.sum())


print(f'A = {A}\nC = {C}\nseeds {seeds}; continental nodes with crust >= 20 km')
print('snapshot | mean age A | C | gap C-A (sd) | seeds C>A | share < 200 Myr A | C')
for tag in tags:
    ra = np.array([read(A, s, tag) for s in seeds]); rc = np.array([read(C, s, tag) for s in seeds])
    gap = rc[:, 0]-ra[:, 0]
    print(f'{tag} | {ra[:, 0].mean():7.1f} | {rc[:, 0].mean():7.1f} | {gap.mean():+6.1f} ({gap.std(ddof=1):5.1f}) | '
          f'{int((gap > 0).sum())} | {ra[:, 1].mean():.3f} | {rc[:, 1].mean():.3f}')
