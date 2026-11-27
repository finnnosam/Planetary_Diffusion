"""Equal-area maps of the same plate arrays and neighbourhood test used in the report."""
from terrain_probe import g,nbr,HERE,LATEST,LONG,REF,ice_mask
from planet_grid import cube
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb

def run():
    x,y=np.meshgrid(np.linspace(-2*np.sqrt(2),2*np.sqrt(2),640),np.linspace(np.sqrt(2),-np.sqrt(2),320))
    th=np.arcsin(np.clip(y/np.sqrt(2),-1,1))
    lo=x*np.pi/(2*np.sqrt(2)*np.maximum(np.cos(th),1e-8))
    la=np.arcsin(np.clip((2*th+np.sin(2*th))/np.pi,-1,1))
    inside=np.abs(lo)<=np.pi
    xyz=np.stack((np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)),-1)
    f,r,c=cube.coordinates(xyz,g.n);side=g.n+1
    index=g.node_of_entry[f*side**2+np.rint(r).astype(int)*side+np.rint(c).astype(int)]
    colors=['#234461','#5793a6','#a6c5bf','#d4bc8c','#538357','#ab7954','#694c55','#d9dce1']
    labs=['Below −3,000 m','−3,000 to −200 m','−200 to 0 m','Land, 0–2,000 m','Low, flat interior','2,000–4,000 m','Above 4,000 m','Reference mask']
    lut=np.array([to_rgb(x) for x in colors])
    inputs=[('Earth • ETOPO',REF/'earth-n217.npz')]+[(f'Seed {d.name[4:]} • 250 Myr',d/'snapshots/t+0000.0.npz') for d in sorted(LATEST.glob('seed*'))]+[('Seed 1001 • 500 Myr',LONG/'seed1001/snapshots/t+0000.0.npz')]
    fig,axes=plt.subplots(4,2,figsize=(10,10.9))
    for i,(ax,(name,path)) in enumerate(zip(axes.flat,inputs)):
        with np.load(path) as s: z=s['elevation_m']
        mask=~ice_mask(g) if i==0 else np.ones(g.size,bool)
        valid=mask.copy(); low=z.copy();hi=z.copy()
        for _ in range(4):
            low=low[nbr].min(1);hi=hi[nbr].max(1);valid=valid[nbr].all(1)
        k=np.select([z<=-3000,z<=-200,z<=0,z>4000,z>=2000],[0,1,2,6,5],default=3)
        k[valid&(low>0)&(z<500)&(hi-low<300)]=4;k[~mask]=7
        rgb=lut[k[index]];rgb[~inside]=1
        ax.imshow(rgb,interpolation='nearest');ax.axis('off');ax.set_title(name,fontsize=11,loc='left',pad=1,color='#203749')
    fig.legend([Patch(color=c) for c in colors],labs,ncol=4,loc='lower center',frameon=False,fontsize=9)
    fig.subplots_adjust(left=.02,right=.98,bottom=.07,top=.98,wspace=.06,hspace=.14)
    fig.savefig(HERE/'terrain-maps.png',dpi=160,facecolor='white')
if __name__=='__main__':
    run()

