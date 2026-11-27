"""Area and coastline contact of small connected features on the n217 graph."""
import numpy as np
from terrain_probe import g,nbr,a,REF,LATEST,HERE,ice_mask
from attribution_probe import components,dump

def run():
    rows=[]
    for name,path in [('Earth',REF/'earth-n217.npz')]+[(p.name,p/'snapshots/t+0000.0.npz') for p in sorted(LATEST.glob('seed*'))]:
        z=np.load(path)['elevation_m'];land=z>0
        keep=~ice_mask(g) if name=='Earth' else np.ones_like(land)
        # Unique axis-edge shore crossings where both ends survive the comparison mask.
        i,j=np.where(g.axis&(np.arange(g.size)[:,None]<nbr));j=nbr[i,j]
        ok=(land[i]!=land[j])&keep[i]&keep[j]
        i,j=i[ok],j[ok]
        r=dict(name=name,shore_crossings=len(i),crossings_per_mkm2=len(i)/(a[keep].sum()/1e12))
        for label,mask in [('land',land),('water',~land)]:
            comps,sizes=components(mask,nbr,a)
            count=np.bincount(comps[mask])
            shore=np.where(mask[i],i,j)
            for n in (1,5,25):
                small=(comps>=0)&(count[np.maximum(comps,0)]<=n)
                r[f'{label}{n}_area']=100*a[small&keep].sum()/a[mask&keep].sum()
                r[f'{label}{n}_cross']=100*np.mean(small[shore])
                r[f'{label}{n}_count']=np.unique(comps[small&keep]).size
        rows.append(r)
    dump('coast-components.csv',rows)
    for k in rows[0]:
        if k!='name':
            print(k,'Earth vs mean (min,max)',round(rows[0][k],3), np.round([np.mean([r[k] for r in rows[1:]]),
                       min(r[k] for r in rows[1:]),max(r[k] for r in rows[1:])],3))
if __name__=='__main__':run()

