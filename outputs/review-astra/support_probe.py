"""Matched-area, common-support and hand-off sensitivity. Repo venv, repository root."""
from terrain_probe import g,a,nbr,HERE,OUT,LATEST,REF,ice_mask,quant
from attribution_probe import components,interiors, fraction,dump
from planet_grid import cube
from planet_earth.crust1 import crust1_cell
from planet_landscape.coupling import PlateHistory
from planet_tectonics.isostasy import sea_level_for_volume
import numpy as np

def run():
    zearth=np.load(REF/'earth-n217.npz')['elevation_m']
    ice=ice_mask(g);good=~ice
    lat=np.degrees(np.arcsin(g.xyz[:,2]));lon=np.degrees(np.arctan2(g.xyz[:,1],g.xyz[:,0]))
    tibet=(lat>25)&(lat<40)&(lon>70)&(lon<105)
    print('Earth high land area >4km inside Himalayan/Tibetan box / globally (%)',fraction(tibet,good&(zearth>4000)))
    target=fraction(zearth>0,good)/100
    cell=crust1_cell(g);w=np.bincount(cell,a,minlength=64800)
    def avg(x):
        return (np.bincount(cell,a*x,minlength=64800)/np.maximum(w,1))[cell]
    Tearth=np.loadtxt(OUT.parent/'data/crust1/crsthk.xyz')[:,2][cell]*1e3
    Searth=np.loadtxt(OUT.parent/'data/crust1/sedthk.xyz')[:,2][cell]*1e3
    obs=np.load(REF/'earth-seafloor-age-n217.npz')
    rows=[]
    for label in ['earth']+[d.name for d in sorted(LATEST.glob('seed*'))]:
        if label=='earth':
            z=zearth;T=Tearth-Searth;dated=(obs['dated_share']>=.99); age=obs['age_myr'];sel=good
        else:
            sn=np.load(LATEST/label/'snapshots/t+0000.0.npz');z=sn['elevation_m'];T=sn['crust_thickness_m']
            age=sn['ocean_age_myr'];dated=~sn['continental'].astype(bool);sel=np.ones(g.size,bool)
        q=quant(z[sel],a[sel],[1-target])[0] if label!='earth' else 0.
        core,lf,_=interiors(z-q,sel)
        lab,sizes=components(z<=0,nbr,a); ocean=lab==np.argmax(sizes)
        dry=(~ocean)&sel
        for _ in range(4):
            dry=dry[nbr].all(1)
        r=dict(run=label,matched_core=core,matched_lowflat=lf,open_sea_core=fraction(dry,sel&(z>0)))
        zone=(avg(T)>=30000)&(avg(z)>=-1000)&(avg(z)<-200)
        r['common_drowned']=fraction(zone&sel,np.ones(g.size,bool))
        for lo,hi in ((20,60),(60,100)):
            m=dated&(age>=lo)&(age<hi)&(z< -1500)&sel
            inside=m[nbr].all(1)&m
            dz=z[nbr].max(1)-z[nbr].min(1)
            r[f'floor_relief_{lo}']=quant(dz[inside],a[inside],[.5])[0]
        if label!='earth':
            # One fine-node location per coarse quad, at its centre (n868 node 4i+2).
            # Use the same sampling code as build_forcing. All cells sampled, area weighted.
            y,x=np.mgrid[0:g.n,0:g.n]+.5
            pt=np.concatenate([cube.directions(f,y,x,g.n).reshape(-1,3) for f in range(6)])
            wf=np.tile(cube.area_weight(y,x,g.n).ravel(),6)
            hist=PlateHistory(LATEST/label/'snapshots')
            h=hist.sample(0.,'elevation_m',pt);og=hist.sample(0.,'orogeny_age_myr',pt)
            v=np.nan_to_num(sn['orogeny_age_myr'],nan=-1.)
            # Each zero corner must be a finite zero to guarantee a wholly active stencil.
            live=cube.sample(g.to_cube((v==0).astype(np.float32))[None],pt)[0]>.99999
            r['active_marked_old_land']=100*wf[(h>0)&live&(og==0)].sum()/wf[h>0].sum()
            r['active_marked_old_high']=100*wf[(h>3000)&live&(og==0)].sum()/wf[h>3000].sum()
        rows.append(r)
    dump('support.csv',rows)
    for k in rows[-1]:
        if k=='run': continue
        print(k,'Earth',rows[0].get(k),'six-seed mean / range',
              np.round([np.mean([r[k] for r in rows[1:]]),min(r[k] for r in rows[1:]),max(r[k] for r in rows[1:])],3))

if __name__=='__main__':
    run()

