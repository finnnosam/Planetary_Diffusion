"""Read the intermediate terms of the shipped return method on a short run.
The observer returns its original arrays unchanged. Run from repository root.
"""
import csv
import argparse
from pathlib import Path
import numpy as np
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.plates import PlateSim

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--steps',type=int,default=10)
    parser.add_argument('--seed',type=int,default=1001)
    parser.add_argument('--out',type=Path,default=HERE/'remap_trace.csv')
    args=parser.parse_args()
    g=CubeSphere(217)
    sim=PlateSim(g,seed=args.seed,params={'remap_rule':'overlap','rift_area':'plates','threads':2})
    original=sim._tent_return
    parts={}
    subduction=sim._subduction_thickening
    erosion=sim._erode
    def record_sub():
        v=subduction()
        parts['sub']=v.copy()
        return v
    def record_erosion(thick):
        v=erosion(thick)
        parts['ero']=v[0].copy()
        return v
    sim._subduction_thickening=record_sub
    sim._erode=record_erosion
    rows=[]
    def observed(k,pl,wn,own,d,name):
        j,out=original(k,pl,wn,own,d,name)
        if name=='thick':
            ws,cs,w=sim._tent[k]
            w=w*g.area[ws,None]
            coverage=np.bincount(cs.ravel(),weights=w.ravel(),minlength=g.size)
            mass=np.bincount(cs.ravel(),weights=(w*d[ws,None]).ravel(),minlength=g.size)
            local=np.where(coverage[j]>0,mass[j]/np.maximum(coverage[j],1e-15),d[wn[j]])
            oc=~pl.cont[j]
            if oc.any():
                corr=out[oc]-local[oc]
                assert np.ptp(corr)<1e-7
                before=pl.thick[j].astype(float)
                floor=sim.p['ocean_thick_m']
                # Previously normal or thicker nodes crossing 6999 m in this call.
                new=oc&(before>=floor-1)&(before+out<floor-1)
                crosses_without=(before+local<floor-1)
                def projected(v):
                    total=np.bincount(cs.ravel(),weights=(w*v[ws,None]).ravel(),minlength=g.size)
                    return np.where(coverage[j]>0,total[j]/np.maximum(coverage[j],1e-15),v[wn[j]])
                sub=projected(parts['sub']) if new.any() else np.zeros(j.size)
                ero=projected(parts['ero']) if new.any() else np.zeros(j.size)
                rows.append(dict(time=sim.t,slot=k,offset_m=corr[0],
                    new_nodes=int(new.sum()),solely_offset=int((new&~crosses_without).sum()),
                    zero_local=int((new&(abs(local)<1e-9)).sum()),
                    formerly_normal=int((new&(abs(before-floor)<1)).sum()),
                    would_survive_without_local_sub=int((new&(before+out-sub>=floor-1)).sum()),
                    would_survive_without_local_ero=int((new&(before+out-ero>=floor-1)).sum()),
                    fallback_nodes=int((new&(coverage[j]<=0)).sum()),
                    loss_km3=float(np.sum(np.maximum(floor-before[new]-out[new],0)*g.area[j[new]])/1e9),
                    oc_area_Mkm2=float(g.area[j[oc]].sum()/1e12)))
        return j,out
    sim._tent_return=observed
    sim.state()
    for step in range(args.steps):
        sim.step()
        st=sim.state()
        print(step+1, float(st['crust_thickness_m'][~st['continental']].min()),flush=True)
        if step%10==9:
            archive=HERE.parents[2]/'outputs/planet-calibration/sweep-20261126-161947/remap_rule-overlap_rift_area-plates'/f'seed{args.seed}/snapshots'/f't{sim.t:+07.1f}.npz'
            with np.load(archive) as saved:
                print('archive match',all(np.array_equal(saved[k],v,equal_nan=True) for k,v in st.items()),flush=True)
    with args.out.open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print({k:sum(r[k] for r in rows) for k in ('new_nodes','solely_offset','zero_local')})


if __name__=='__main__':
    main()
