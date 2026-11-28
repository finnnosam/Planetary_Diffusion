"""Separate inland water from the open ocean. Check plateau support on 1-degree cells."""
import csv
from pathlib import Path
import numpy as np
from spatial import A,G,CAL,HERE,ROOT,REF,quant
from planet_earth.reference import ice_mask
from planet_earth.crust1 import OCEANIC,crust1_cell
from planet_tectonics.graphops import components,nearest_label

cell=crust1_cell(G)
weights=np.bincount(cell,weights=A,minlength=64800)


def assess(name,z,T,cont,domain):
    r=dict(name=name)
    wet=z<=0
    land=~wet&domain
    labels=components(wet,G.nbr)
    seas=np.bincount(labels[wet],weights=A[wet])
    salt=labels==seas.argmax()
    for kind,mask in [('wet',wet),('open',salt)]:
        _,dist=nearest_label(np.where(mask,0,-1),G.nbr,G.dist)
        for km in (300,500,1000):
            far=land&(dist>=km*1e3)
            r[f'{kind}_far{km}_pct']=100*A[far].sum()/A[land].sum()
            r[f'{kind}_low_among_far{km}_pct']=100*A[far&(z<200)].sum()/A[far].sum()
    # Both variables averaged to the same 1-degree cells; >50% continental coverage.
    def avg(v): return np.bincount(cell,weights=A*v,minlength=64800)/np.maximum(weights,1.)
    t=avg(T)
    h=avg(z)
    m=(avg(cont.astype(float))>.5)&(avg(domain.astype(float))>.5)&(weights>0)
    for lo,hi in [(45,50),(50,55),(55,60),(60,65),(65,70),(70,80)]:
        valid=m&(t>=lo*1000)&(t<hi*1000)
        r[f'area_{lo}_{hi}_Mkm2']=weights[valid].sum()/1e12
        r[f'z_{lo}_{hi}_m']=quant(h[valid],weights[valid]) if valid.any() else np.nan
    return r


def main():
    z=np.load(REF/'earth-n217.npz')['elevation_m']
    t=np.loadtxt(ROOT/'data/crust1/crsthk.xyz')[:,2][cell]*1000
    typ=np.loadtxt(ROOT/'data/crust1/CNtype1-1.txt',dtype=str).ravel()[cell]
    rows=[assess('Earth',z,t,~np.isin(typ,list(OCEANIC)),~ice_mask(G))]
    for short,arm in [('A','remap_rule-nearest_rift_area-world'),('C','remap_rule-overlap_rift_area-plates')]:
        for p in sorted((CAL/'sweep-20261126-161947'/arm).glob('seed*/snapshots/t+0000.0.npz')):
            with np.load(p) as f:
                rows.append(assess(short+'/'+p.parents[1].name,f['elevation_m'],f['crust_thickness_m']+
                                   f['sediment_bulk_m'],f['continental'],np.ones(G.size,bool)))
    with (HERE/'interiors.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    for key in list(rows[0])[1:]:
        summary=[rows[0][key]]
        for arm in (rows[1:7],rows[7:]):
            finite=[r[key] for r in arm if np.isfinite(r[key])]
            summary.append(np.median(finite) if finite else np.nan)
        print(key,*(round(x,2) for x in summary))


if __name__=='__main__':
    main()
