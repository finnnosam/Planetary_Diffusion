"""Exactly sum the archived forcing with landscape surface processes switched off.

No simulation or parameter changes. The flexure operator is linear. Use the normal
fine-grid geometry without allocating its neighbour graph.
"""
from pathlib import Path
import csv

import numpy as np
from planet_grid import cube
from planet_grid.sphere_graph import CubeSphere, _quad_areas, EARTH_RADIUS_M
from planet_grid.sphere_flexure import SphereFlexure
from planet_landscape.graph_landscape import DEFAULTS
from planet_earth.reference import ice_mask
from planet_tectonics.graphops import components

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CHAIN=ROOT/'outputs/review-astra/runs/Chain-360f894f'
REPLAY=ROOT/'outputs/planet-landscape/age-replay-20261123'


def geometry(n):
    side=n+1
    members,ids,owners,_=cube.edge_groups(n)
    canon=np.arange(6*side**2)
    canon[members]=owners[ids]
    unique=np.unique(canon)
    node=np.searchsorted(unique,canon).reshape(6,side,side)
    xyz=np.empty((len(unique),3))
    area=np.zeros(len(unique))
    yy,xx=np.indices((side,side))
    for f,idx in enumerate(node):
        pos=cube.directions(f,yy,xx,n)
        xyz[idx]=pos
        q=_quad_areas(pos)*EARTH_RADIUS_M**2/4
        for y,x in [(0,0),(0,1),(1,0),(1,1)]:
            area[idx[y:y+n,x:x+n]]+=q
    assert np.isclose(area.sum(),4*np.pi*EARTH_RADIUS_M**2)
    return xyz,area


def draw(fields,coarse_grid):
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    c=coarse_grid
    earth=np.load(ROOT/'outputs/planet-reference/earth-n217.npz')['elevation_m']
    parent=np.load(CHAIN/'plates/seed1001/snapshots/t+0000.0.npz')['elevation_m']
    def centre(z,valid):
        m=(z>3000)&valid
        lab=components(m,c.nbr)
        winner=np.bincount(lab[m],weights=c.area[m]).argmax()
        ids=lab==winner
        v=np.average(c.xyz[ids],axis=0,weights=c.area[ids])
        return v/np.linalg.norm(v)
    centers=[centre(earth,~ice_mask(c)),centre(parent,np.ones(c.size,bool))]
    n=868
    members,ids,owners,_=cube.edge_groups(n)
    canon=np.arange(6*(n+1)**2)
    canon[members]=owners[ids]
    node=np.searchsorted(np.unique(canon),canon).reshape(6,n+1,n+1)
    def window(v,center):
        east=np.cross([0,0,1.],center);east/=np.linalg.norm(east)
        north=np.cross(center,east)
        xx,yy=np.meshgrid(np.linspace(-1400,1400,351),np.linspace(-950,950,239))
        r=np.hypot(xx,yy)
        u=r/(EARTH_RADIUS_M/1000)
        p=np.cos(u)[...,None]*center+(xx[...,None]*east+yy[...,None]*north)*(np.sinc(u/np.pi)/(EARTH_RADIUS_M/1000))[...,None]
        face,y,x=cube.coordinates(p,n)
        return v[node[face,np.rint(y).astype(int),np.rint(x).astype(int)]]
    fig,axs=plt.subplots(2,2,figsize=(9.2,6.8))
    cmap=ListedColormap(['#a6c8d8','#c6dcbf','#dce5b3','#efe4a5','#dcbd78','#be925a','#925e3e','#674135','#382b32'])
    levels=[-12000,0,200,500,1000,2000,3000,4000,5000,12000]
    for ax,(field,label,idx) in zip(axs.flat,[('earth','Earth · Tibetan Plateau',0),('parent','Planet · plate surface',1),
            ('forcing_only','Same planet · forcing alone',1),('plate','Same planet · finished landscape',1)]):
        a=window(fields[field],centers[idx])
        im=ax.imshow(a,origin='lower',extent=(-1400,1400,-950,950),cmap=cmap,norm=BoundaryNorm(levels,cmap.N),interpolation='nearest')
        ax.set_title(label,loc='left',fontsize=12)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
    fig.subplots_adjust(left=.025,right=.98,top=.95,bottom=.14,hspace=.14,wspace=.07)
    cb=fig.colorbar(im,cax=fig.add_axes([.09,.073,.82,.025]),orientation='horizontal',
                    ticks=levels[1:-1],spacing='uniform')
    cb.set_label('Elevation relative to sea level (m)    ·    Each panel is 2,800 × 1,900 km')
    fig.savefig(HERE/'high-country.png',dpi=145,bbox_inches='tight')
    plt.close(fig)


def main():
    xyz,a=geometry(868)
    f=CHAIN/'landscape/forcing'
    z0=np.load(f/'z0_m.npy')+np.random.default_rng(0).normal(0,1.,len(a))
    load=np.load(f/'load_m.npy',mmap_mode='r').sum(axis=0,dtype='float64')
    surface=np.load(f/'surface_m.npy',mmap_mode='r').sum(axis=0,dtype='float64')
    flex=SphereFlexure(xyz,a,EARTH_RADIUS_M,DEFAULTS['Te'])
    noerosion=z0+load+surface-flex.deflection(load)
    c=CubeSphere(217)
    with np.load(CHAIN/'plates/seed1001/snapshots/t+0000.0.npz') as snap:
        coarse=c.to_cube(snap['elevation_m'])[None]
    endpoint=np.concatenate([cube.sample(coarse,points)[0] for points in np.array_split(xyz,32)])
    earth=np.load(ROOT/'outputs/planet-reference/earth-n868.npz')['elevation_m']
    exclude=ice_mask(type('Coords',(),{'xyz':xyz})())
    fields={'earth':earth,'parent':endpoint,'initial':z0,'forcing_only':noerosion}
    for arm in ['legacy','plate']:
        fields[arm]=np.load(REPLAY/arm/'landscape/present.npz')['elevation_m']
    rows=[]
    for tag,z in fields.items():
        land=(z>0)&(~exclude if tag=='earth' else True)
        ind=np.flatnonzero(land); ix=ind[np.argsort(z[ind])]
        cw=np.cumsum(a[ix]); q=z[ix[np.searchsorted(cw,np.array([.5,.9,.99])*cw[-1])]]
        row=dict(tag=tag,p50=q[0],p90=q[1],p99=q[2],pct_land_gt3km=100*a[land&(z>3000)].sum()/a[land].sum())
        if tag!='earth':
            for low,high in [(0,500),(500,2000),(2000,3000),(3000,np.inf)]:
                m=(endpoint>low)&(endpoint<=high)
                row[f'mean_on_parent_{low}m']=np.average(z[m],weights=a[m])
        rows.append(row)
    with (HERE/'replay.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[-1]); w.writeheader(); w.writerows(rows)
    for row in rows: print(row)
    draw(fields,c)


if __name__=='__main__':
    main()
