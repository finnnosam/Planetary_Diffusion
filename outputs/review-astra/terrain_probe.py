"""Run with the repository venv from any directory. Outputs are beside this script."""
from pathlib import Path
import csv
import json
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_earth.reference import ice_mask

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
CAL = OUT / 'planet-calibration'
REF = OUT / 'planet-reference'
LATEST = CAL / 'sweep-20261122-150810/combo-01'
PREV = CAL / 'sweep-20261121-220857'
LONG = OUT / 'planet-plates/long-500-seed1001-20261122'
g = CubeSphere(217)
a = g.area
valid = np.ones(g.size, bool)
# Draw uniform-area comparisons by weighting each centre. Fixed coarse support:
# extrema over the centre and its k-step 8-neighbour graph ball.
nbr = g.nbr.copy()
bad = nbr < 0
nbr[bad] = np.broadcast_to(np.arange(g.size)[:, None], nbr.shape)[bad]

def quant(z, w, qs):
    order = np.argsort(z)
    cw = np.cumsum(w[order])
    return z[order][np.searchsorted(cw, np.asarray(qs)*cw[-1])]

def describe(z, mask, full):
    l = (z > 0) & mask
    o = (z <= 0) & mask
    def pct(m, base):
        return 100*a[m & base].sum()/a[base].sum()
    def mean(x, m):
        return np.average(x[m], weights=a[m])
    row = dict(land_pct=pct(l, mask), mean_land=mean(z,l), mean_sea=mean(z,o),
               low200=pct(z<200,l), mid500_2000=pct((z>=500)&(z<2000),l),
               high3000=pct(z>3000,l), high4000=pct(z>4000,l),
               shelf200=pct(z>=-200,o), abyss5000=pct(z< -5000,o))
    for key,m in [('land',l),('sea',o)]:
        qs=[.1,.5,.9,.99]
        row.update(zip([f'{key}{int(q*100)}' for q in qs],quant(z[m],a[m],qs)))
    if full:
        low,high=z.copy(),z.copy()
        ok=mask.copy()
        for k in range(1,5):
            low=low[nbr].min(1)
            high=high[nbr].max(1)
            ok=ok[nbr].all(1)
            if k in (1,4):
                dry=ok & (low>0)
                wet=ok & (high<=0)
                relief=high-low
                row[f'landinterior{k}']=pct(dry,l)
                row[f'flat{k}']=pct(relief<300,dry)
                row[f'lowflat{k}']=pct((relief<300)&(z<500),dry)
                row[f'rugged{k}']=pct(relief>1500,dry)
                for tag,sel in [('land',dry),('sea',wet)]:
                    row[f'{tag}relief{k}']=quant(relief[sel],a[sel],[.5])[0]
    return row

def run():
    rows=[]
    z=np.load(REF/'earth-n217.npz')['elevation_m']
    for tag,mask in [('earth-noice',~ice_mask(g)),('earth-all',valid)]:
        rows.append(dict(run=tag,seed=0,time=0,**describe(z,mask,True)))
    for tag,path,times in [('repeat',LATEST,[-40,-30,-20,-10,0]),('oldOn',PREV/'rift_traction-True',[0]),('oldOff',PREV/'rift_traction-False',[0]),('long',LONG,[-500,-450,-400,-350,-300,-250,-200,-150,-100,-50,-40,-30,-20,-10,0])]:
        for d in sorted(path.glob('seed*')):
            stats=json.loads((d/'stats.json').read_text())
            for t in times:
                with np.load(d/'snapshots'/f't{t:+07.1f}.npz') as sn:
                    assert sn['time_myr']==t
                    z=sn['elevation_m']
                    r=describe(z,valid, t==0 or tag=='long' and t%100==0)
                    c=sn['continental'].astype(bool)
                    age=sn['ocean_age_myr']
                    r['floor_age']=np.average(age[~c],weights=a[~c])
                    r['floor_under20']=100*a[(~c)&(age<20)].sum()/a[~c].sum()
                    r['floor_prod10']=a[(~c)&(age<10)].sum()/1e12/10
                    r['thick55']=100*a[c & (sn['crust_thickness_m']+sn['sediment_bulk_m']>55000)].sum()/a[c].sum()
                    r['sed_cont_km']=np.average(sn['sediment_bulk_m'][c],weights=a[c])/1000
                    r['sed_ocean_km']=np.average(sn['sediment_bulk_m'][~c],weights=a[~c])/1000
                    q=next(item for item in stats['series'] if item['time_myr']==t)
                    assert abs(r['mean_land']-q['land_mean_elevation_m'])<1.
                    for k in ['plates','plates_over_1mkm2','continental_fraction','continental_volume_km3','continental_sediment_km3','sea_level_on_datum_m']:
                        r[k]=q.get(k)
                    rows.append(dict(run=tag,seed=d.name,time=t,**r))
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (HERE/'terrain.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

if __name__=='__main__':
    run()

