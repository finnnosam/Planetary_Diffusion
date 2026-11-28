"""Read plateau capacity on every archived 10 Myr snapshot (no new realization)."""
import csv
import numpy as np
from spatial import A,G,CAL,HERE,REF
from planet_earth.reference import ice_mask
from planet_tectonics.graphops import components


def measure(z,valid):
    m=(z>4000)&valid
    if not m.any():
        return 0.,0.
    lab=components(m,G.nbr)
    return np.bincount(lab[m],weights=A[m]).max()/1e12, A[m].sum()/1e12


def main():
    rows=[]
    e=measure(np.load(REF/'earth-n217.npz')['elevation_m'],~ice_mask(G))
    for tag,dirname in [('A','remap_rule-nearest_rift_area-world'),('C','remap_rule-overlap_rift_area-plates')]:
        for seed in sorted((CAL/'sweep-20261126-161947'/dirname).glob('seed*')):
            by_seed=[]
            for t in range(-250,1,10):
                with np.load(seed/'snapshots'/f't{t:+07.1f}.npz') as f:
                    largest,total=measure(f['elevation_m'],np.ones(G.size,bool))
                by_seed.append(dict(arm=tag,seed=seed.name,time=t,largest_Mkm2=largest,total_Mkm2=total))
            rows+=by_seed
            print(tag,seed.name, max(by_seed,key=lambda r:r['largest_Mkm2']))
    print('Earth',e)
    with (HERE/'plateau_history.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)


if __name__=='__main__':
    main()
