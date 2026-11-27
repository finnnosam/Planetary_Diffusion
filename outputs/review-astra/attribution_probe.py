"""Post-first-read checks. Run at repository root with its Python venv."""
from terrain_probe import (OUT, HERE, LATEST, PREV, REF, LONG, quant, g, a,
                           nbr, valid, ice_mask)
from planet_earth.crust1 import crust1_cell
from planet_tectonics import mantle_flow as mf
from planet_tectonics.isostasy import sea_level_for_volume
from planet_grid import cube
from pathlib import Path
import h5py
import numpy as np
import json, csv
from numba import njit

@njit
def components(m, adj, area):
    lab=np.full(m.size,-1)
    queue=np.empty(m.size,np.int32)
    sizes=np.zeros(m.size)
    n=0
    for i in range(m.size):
        if not m[i] or lab[i]>=0: continue
        queue[0]=i; lab[i]=n; end=1; pos=0
        while pos<end:
            k=queue[pos];pos+=1;sizes[n]+=area[k]
            for v in adj[k]:
                if m[v] and lab[v]<0:
                    lab[v]=n;queue[end]=v;end+=1
        n+=1
    return lab,sizes[:n]

def fraction(m,base):
    return 100*a[m&base].sum()/a[base].sum()

def interiors(z,mask):
    low=z.copy(); high=z.copy(); ok=mask.copy()
    for _ in range(4):
        low=low[nbr].min(1); high=high[nbr].max(1); ok=ok[nbr].all(1)
    dry=(low>0)&ok
    return [fraction(dry,mask&(z>0)),fraction(dry&(z<500)&(high-low<300),mask&(z>0)),
            quant((high-low)[ok&(high< -1500)],a[ok&(high< -1500)],[.5])[0]]
def geog(z,mask):
    # Connectivity on the whole native grid, reports use the same local mask as terrain.csv.
    lab,sizes=components(z<=0,nbr,a)
    sea=lab==np.argmax(sizes)
    _,mountains=components((z>4000)&mask,nbr,a)
    return dict(zip(('dryCore','lowFlatCore','deep_relief'),interiors(z,mask)),
                disconnected_wet=fraction(~sea,mask&(z<=0)),
                max_high_area=mountains.max()/1e12 if mountains.size else 0)
def clustering(T):
    c=T>=20000
    lm=10
    alm=mf.point_alm(g.xyz,(c-np.average(c,weights=a))*a/g.radius**2,lm)
    p=mf.spectrum(alm,lm)[1:]
    return p[0]/p.sum(),np.linalg.norm((a[c,None]*g.xyz[c]).sum(0))/a[c].sum()

def dump(name,rows):
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader();w.writerows(rows)

def observed_rows():
    root=OUT.parent; data=root/'data'; cell=crust1_cell(g)
    mask=~ice_mask(g); z=np.load(REF/'earth-n217.npz')['elevation_m']
    T,S=[np.loadtxt(data/'crust1'/f)[:,2][cell]*1000 for f in ('crsthk.xyz','sedthk.xyz')]
    print('Earth continent (proper total vs total + cover), area %, degree1, resultant')
    for v in (T,T+S):
        print(fraction(v>=20000,valid),clustering(v))
    print('Earth raw and common-support crust thresholds. Share of nonice / whole globe')
    def lump(x):
        w=np.bincount(cell,a,minlength=64800)
        return (np.bincount(cell,a*x,minlength=64800)/np.maximum(w,1))[cell]
    for height in (z,lump(z)):
        for threshold in (27000,30000,33000):
            m=(T-S>=threshold)&(height>=-1000)&(height< -200)
            print(threshold, fraction(m,mask), fraction(m&mask,valid))
    grd=next(data.rglob('degree_1_to_30.grd'))
    def grid_h(lon,lat):
        with h5py.File(grd) as f:
            return f['z'][:][np.rint(lat+90).astype(int),np.rint(lon%360).astype(int)]*1000
    lat=np.rad2deg(np.arcsin(g.xyz[:,2]));lon=np.rad2deg(np.arctan2(g.xyz[:,1],g.xyz[:,0]))
    h=grid_h(lon,lat)
    for typ in ['continental_measurements_seismic','oceanic_measurements_circles']:
        obs=np.loadtxt(next(grd.parent.glob(typ+'*')))
        x=grid_h(*obs[:,:2].T);y=obs[:,2]*1000
        print('Holdt raw points vs its own grid (descriptive in-sample)',typ,len(y),
              'r',np.corrcoef(x,y)[0,1],'RMS difference m',np.mean((x-y)**2)**.5)
    m=(T-S>=30000)&(z>=-1000)&(z< -200)&mask
    pts=np.loadtxt(next(grd.parent.glob('continental_measurements_seismic*')))
    la,lo=np.deg2rad(pts[:,1]),np.deg2rad(pts[:,0])
    xyz=np.stack((np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)),1)
    f,yy,xx=cube.coordinates(xyz,g.n); side=g.n+1
    idx=g.node_of_entry[f*side**2+np.rint(yy).astype(int)*side+np.rint(xx).astype(int)]
    print('drowned normal Earth direct points n, 1-degree cells, global area %, Holdt grid mean m',
          m[idx].sum(),np.unique(cell[m]).size,fraction(m,valid),np.average(h[m],weights=a[m]))
    return dict(run='Earth',seed=0,time=0,**geog(z,mask))

def run():
    rows=[observed_rows()];null=[]
    # Rigid cube symmetries preserve node areas and the spectrum exactly. Subtract saved
    # realised loaded field, rotate it, restore volume. This is a diagnostic, not rerunning feedback.
    from itertools import permutations,product
    idx=np.arange(g.size)
    rots=[]
    for perm in permutations(range(3)):
        for sign in product([-1,1],repeat=3):
            mat=np.eye(3)[:,perm]*sign
            if np.linalg.det(mat)<0 or np.array_equal(mat,np.eye(3)): continue
            f,y,x=cube.coordinates(g.xyz@mat,217); side=g.n+1
            j=g.node_of_entry[f*side**2+np.rint(y).astype(int)*side+np.rint(x).astype(int)]
            rots.append(j)
    for tag,arm in [('off',PREV/'rift_traction-False'),('on',PREV/'rift_traction-True'),('repeat',LATEST)]:
        for sd in sorted(arm.glob('seed*')):
            for t in (-40,-30,-20,-10,0):
                with np.load(sd/'snapshots'/f't{t:+07.1f}.npz') as snap:
                    z=snap['elevation_m']; zd=z+float(snap['sea_level_m'])
                    T=snap['crust_thickness_m']; total=T+snap['sediment_bulk_m']; flag=snap['continental'].astype(bool)
                    typ=flag&(total>=30000); normal=(T>=30000)
                    dr=normal&(z>=-1000)&(z< -200); dyn=snap['dynamic_topography_m']
                    l1,res=clustering(total)
                    row=dict(run=tag,seed=sd.name,time=t,L1=l1,resultant=res,
                             datum10=quant(zd[typ],a[typ],[.1])[0],datum50=quant(zd[typ],a[typ],[.5])[0],
                             drownedDatum=fraction(zd<0,typ),drownedNormal=fraction(dr,valid),
                             dynDrowned=np.average(dyn[dr],weights=a[dr]))
                    if t==0:
                        row.update(geog(z,valid))
                    if t==0 and tag=='repeat':
                        vol=np.sum(a*np.maximum(-z,0))
                        bare=zd-dyn
                        for j,field in enumerate([np.zeros_like(dyn),dyn]+[dyn[k] for k in rots]):
                            b=bare+field
                            sea=sea_level_for_volume(b,a,vol)
                            v=b-sea;m=normal&(v>=-1000)&(v< -200)
                            null.append(dict(seed=sd.name,rotation=j-1,drowned=fraction(m,valid),field=np.average(field[m],weights=a[m]),
                                             land=fraction(v>0,valid),sea=sea))
                        v=bare-sea_level_for_volume(bare,a,vol)
                        row['noFieldCore'],row['noFieldLowFlat'],_=interiors(v,valid)
                    rows.append(row)
    dump('attribution.csv',rows);dump('rotations.csv',null)

if __name__=='__main__':
    run()

