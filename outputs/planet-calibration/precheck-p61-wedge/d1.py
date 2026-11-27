import numpy as np, json
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_reviews.crust_review import column_thickness_m, run_params
from planet_tectonics.plates import PARAMS
g=CubeSphere(217); a=g.area
for arm in ('control','wedge'):
    d=f'/home/claude/pre/{arm}/seed1010'; p=run_params(d)
    for t in ('-0040.0','+0000.0'):
        s=np.load(f'{d}/snapshots/t{t}.npz'); z=s['elevation_m'].astype(float); bt=s['boundary_type']
        T=column_thickness_m(s,p); edge=T>=20e3
        dc=gl.distance_to(edge,g.nbr,g.dist,1e6); far=gl.distance_to(bt!=0,g.nbr,g.dist,3e5)>3e5
        km,dep=np.array(PARAMS['margin_wedge_profile']).T
        top=np.where(~edge&far&(dc<=1e6),np.interp(dc/1e3,km,dep),-1e9)
        on=top>-1e8
        ocean=gl.connected_ocean(z,z<-1500.,g.nbr)
        dl=gl.distance_to(~ocean,g.nbr,g.dist,150e3)
        beyond=ocean&(dl>150e3)
        near_edge=~edge&(dc<=300e3)&ocean
        print(arm,t,'wedge-eligible %.1f Mkm2, room %.1f Mkm3; ocean beyond shelf within 300 km of the edge %.1f Mkm2, of it eligible %.1f; of it within 300 km of a boundary %.1f'%(
            a[on].sum()/1e12, (np.clip(top-z,0,None)*a)[on].sum()/1e15, a[near_edge&beyond].sum()/1e12, a[near_edge&beyond&on].sum()/1e12, a[near_edge&beyond&~far].sum()/1e12))
        # where is the shelf zone edge relative to the continental edge
        print('   coast-to-edge: shelf-zone (<=150 km of land) share of ocean within 0-50 km of the edge: %.2f'%(a[ocean&(dc>0)&(dc<=50e3)&(dl<=150e3)].sum()/max(a[ocean&(dc>0)&(dc<=50e3)].sum(),1)))
        # routing: overflow ends at local lows; where are the lows?
        sea_in=np.zeros(g.size); coast=ocean&(dl<=50e3)&np.isin(np.arange(g.size), np.unique(g.nbr[~ocean][g.nbr[~ocean]>=0]))
        sea_in[coast]=1e12
        dp,left=gl.marine_sediment_kept(z,ocean,sea_in,g.nbr,g.dist,a,10.,1.2e-3,dl,150e3,top)
        wedge=(np.where(dl>150e3,dp,0)*a).sum()
        L=left>0
        print('   test supply %.0f Mkm3: on the shelf zone %.0f, wedge %.0f, left %.0f; left nodes: %d, %.0f%% within 300 km of a boundary, mean depth %.0f m, eligible %.0f%%'%(
            sea_in.sum()/1e15, (np.where(dl<=150e3,dp,0)*a).sum()/1e15, wedge/1e15, left.sum()/1e15, L.sum(), 100*(left[L&~far].sum()/left.sum()), np.average(z[L],weights=left[L]), 100*left[L&on].sum()/left.sum()))
        z2,wv=gl.fill_and_spill(z+dp,left,g.nbr,a,top,return_wedge=True)
        print('   spill path: of the left %.2f Mkm3, on the wedge %.2f'%(left.sum()/1e15, wv.sum()/1e15))
        big=left*50
        z2,wv=gl.fill_and_spill(z+dp,big,g.nbr,a,top,return_wedge=True)
        print('   x50 left %.1f Mkm3, on the wedge %.1f'%(big.sum()/1e15, wv.sum()/1e15))
