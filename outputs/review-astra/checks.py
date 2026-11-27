"""Reproduces ledger changes, repeated-arm uncertainty and fixed-column checks."""
import csv,json
import numpy as np
from terrain_probe import g,a,nbr,HERE,REF,OUT,PREV,LONG,LATEST,quant,ice_mask
from planet_tectonics.isostasy import sea_level_for_volume
from planet_tectonics import mantle
from planet_landscape import coupling
from planet_grid import cube
from planet_grid.fields import nearest_nodes
from planet_earth.crust1 import crust1_cell
from itertools import permutations, product
from pathlib import Path
import rasterio

def load_csv(name):
    with (HERE/name).open() as f:
        return list(csv.DictReader(f))

def reduction():
    rows=load_csv('attribution.csv')
    rng=np.random.default_rng(1139)
    draw=rng.integers(6,size=(20000,6))
    for key in ('L1','datum10','datum50','drownedNormal'):
        for window in (False,True):
            vs=[]
            for run in ('off','on','repeat'):
                vals=[np.mean([float(x[key]) for x in rows if x['seed']==f'seed{s}' and x['run']==run and
                               (window or x['time']=='0')]) for s in range(1001,1007)]
                vs.append(np.array(vals))
            d=vs[1]-vs[0]; delta=vs[2]-vs[0]
            print(key,'window',window,'means',np.round(np.mean(vs,1),3),
                  'adoption minus off, paired bootstrap 95%',np.round(np.percentile(d[draw].mean(1),[2.5,97.5]),3),
                  'repeat minus off',np.round(np.percentile(delta[draw].mean(1),[2.5,97.5]),3))
    print('old vs new effective parameter differences')
    prev=json.loads((PREV/'rift_traction-True/seed1001/stats.json').read_text())['params_effective']
    now=json.loads((LATEST/'seed1001/stats.json').read_text())['params_effective']
    print([(k,prev.get(k),now.get(k)) for k in sorted(prev|now) if prev.get(k)!=now.get(k)])
    print('long run by halves. Rates km3/yr, area contributions percentage points')
    data=json.loads((LONG/'seed1001/stats.json').read_text())
    S={x['time_myr']:x for x in data['series']}
    for start,end in [(-500,-250),(-250,0),(-100,0)]:
        p,q=S[start],S[end]
        dt=(end-start)*1e6
        er=lambda k:(q['erosion_totals_km3'].get(k,0)-p['erosion_totals_km3'].get(k,0))/dt
        net=(q['continental_volume_km3']-p['continental_volume_km3'])/dt
        terms={k:round(er(k),3) for k in ['continental_crust_eroded_m3','subduction_on_continents_m3','sediment_accreted_m3','collision_sediment_to_crust_m3','arc_accreted_on_continent_m3','eclogite_foundered_m3','volcanic_added_m3']}
        print(start,end,'crust',round(net,3),'sed',(q['continental_sediment_km3']-p['continental_sediment_km3'])/dt,terms)
        print({k:round(100*(v-p['area_ledger_frac'].get(k,0)),3) for k,v in q['area_ledger_frac'].items()})
        print('events',{k:v-p['events'].get(k,0) for k,v in q['events'].items()})
    for tt in (-500,-250,0):
        sn=np.load(LONG/'seed1001/snapshots'/f't{tt:+07.1f}.npz')
        print('slab record, mass by layer sum kg',tt,sn['slab_record_kg'],sn['slab_mass_by_layer_kg'].sum())
    for base,tag in [(LONG,-500),(LATEST,-250)]:
        z=np.load(base/'seed1001/snapshots'/f't{tag:+07.1f}.npz')
        print('seed1001 start',tag,np.average(z['elevation_m'][z['elevation_m']>0],weights=a[z['elevation_m']>0]))
    rows=load_csv('terrain.csv')
    for k in ['land_pct','floor_prod10','high4000','mean_land']:
        print(k,'repeat per epoch',
              [round(np.mean([float(r[k]) for r in rows if r['run']=='repeat' and r['time']==str(t)]),3)
                   for t in [-40,-30,-20,-10,0]])

def support():
    valid=~ice_mask(g);z=np.load(REF/'earth-n217.npz')['elevation_m']
    V=np.loadtxt(OUT.parent/'data/mantle/vd2020_sph_coeff_high_accuracy_spot.dat')
    scale=2300/3300*1e3/np.sqrt(4*np.pi)
    print('VD mean map / posterior rms m, posterior variance fraction',
          scale*np.linalg.norm(V[:,2]),scale*np.linalg.norm(V[:,2:4]),np.sum(V[:,3]**2)/np.sum(V[:,2:4]**2))
    # Pixel-to-cube independent spot re-aggregation. Read full narrow strips only around each site.
    with rasterio.open(OUT.parent/'data/etopo_2022_30s_surface.tif') as ds:
        for lon,lat in [(0,-80),(-42,74),(90,30),(20,0),(-70,-22)]:
            rr,cc=ds.index(lon,lat)
            r=np.arange(rr-150,rr+151);c=np.arange(cc-560,cc+561)
            la,lo=np.meshgrid(np.deg2rad(90-(r+.5)/120),np.deg2rad(-180+(c+.5)/120),indexing='ij')
            vec=np.stack((np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)),-1)
            n=nearest_nodes(g,vec.reshape(-1,3)).reshape(la.shape);node=n[len(r)//2,len(c)//2];m=n==node
            v=ds.read(1,window=((r[0],r[-1]+1),(c[0],c[-1]+1)))
            print('raster',lon,lat,'point / pixel mean / cache',round(float(v[len(r)//2,len(c)//2]),1),
                  round(np.average(v[m],weights=np.cos(la[m])),2),round(float(z[node]),2))

def population():
    rates=[]; sets=[]
    for d in sorted(LATEST.glob('seed*')):
        sn=np.load(d/'snapshots/t+0000.0.npz')
        T=sn['crust_thickness_m'];z=sn['elevation_m'];dyn=sn['dynamic_topography_m']
        sea=sn['sea_level_m'];base=z+sea-dyn; vol=np.sum(a*np.maximum(-z,0))
        bare=base-sea_level_for_volume(base,a,vol)
        m=T>=30000; select=lambda x:m&(x>=-1000)&(x< -200)
        orig=select(z); new=select(bare)
        rates.append([100*a[orig&~new].sum()/a[orig].sum(),100*a[new&~orig].sum()/a[new].sum(),
                      *[100*a[m&(v<cut)].sum()/a[m].sum() for v in (z,bare) for cut in (-1000,-200,0)],
                      *[100*a[(T>=h)&(z>=-1000)&(z< -200)].sum()/a.sum() for h in (27000,30000,33000)]])
        og=sn['orogeny_age_myr'];land=z>0
        # Coarse-node age and all zero-age interiors. An all-finite, all-zero chart patch
        # remains exactly zero even after bilinear sampling onto the fine grid.
        zero=(og==0); zi=zero[nbr].all(1)&zero
        sets.append([100*a[zero&land].sum()/a[land].sum(),100*a[zi&land].sum()/a[land].sum(),
                     100*a[zero&(z>3000)].sum()/a[z>3000].sum()])
    print('band departures / entrants %, normal crust below -1000,-200,0 before & after %, threshold sweep %',np.round(np.mean(rates,0),3))
    print('zero-age land %, interiors all zero-age neighbours %, zero-age on >3km land %',np.round(np.mean(sets,0),3))
    print('fine-grid response for a uniform observed age=0 patch',coupling.erodibility(np.where(np.array([0.])>0,0.,1e4)))

if __name__=='__main__':
    reduction();support();population()

