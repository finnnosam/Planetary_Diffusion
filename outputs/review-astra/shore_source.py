"""What makes a single-cell island or hole stand apart from its immediate neighbours?"""
from terrain_probe import g,nbr,a,REF,LATEST,HERE,ice_mask,quant
from attribution_probe import dump
import numpy as np

def run():
    rows=[]
    for name,p in [('Earth',REF/'earth-n217.npz')]+[(d.name,d/'snapshots/t+0000.0.npz') for d in sorted(LATEST.glob('seed*'))]:
        sn=np.load(p); z=sn['elevation_m']; land=z>0
        fields={'z':z}
        if name!='Earth':
            fields.update(dyn=sn['dynamic_topography_m'],air=sn['dynamic_topography_air_m'],crust=sn['crust_thickness_m'],sed=sn['sediment_bulk_m'])
        mask=~ice_mask(g) if name=='Earth' else np.ones(g.size,bool)
        for tag,v in [('island',land),('hole',~land)]:
            m=v&mask&(~v[nbr] | (g.nbr<0)).all(1)&mask[nbr].all(1)
            r=dict(name=name,feature=tag,n=m.sum())
            for k,f in fields.items():
                delta=f[m]-np.nanmedian(np.where(g.nbr[m]>=0,f[nbr[m]],np.nan),axis=1)
                r[k]=quant(delta,a[m],[.5])[0]
            if name!='Earth':
                c=sn['continental'].astype(bool)
                r['inside_continent']=np.average(c[m]&c[nbr[m]].all(1),weights=a[m])
                nodyn=z-sn['dynamic_topography_m']
                r['survives_fixed_sea']=np.average(((nodyn[m]>0)==land[m]) & ((nodyn[nbr[m]]>0)==land[nbr[m]]).all(1),weights=a[m])
            rows.append(r)
    dump('shore-source.csv',rows)
    for tag in ('island','hole'):
        rs=[r for r in rows if r['feature']==tag and r['name']!='Earth']
        print(tag, {k:round(np.mean([r[k] for r in rs]),3) for k in rs[0] if k not in ('name','feature')},'Earth', next(r for r in rows if r['feature']==tag and r['name']=='Earth'))
if __name__=='__main__':run()

