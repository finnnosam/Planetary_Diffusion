"""Use the archived annual forcing path and basin observations, with no new model run."""
import csv
import json
from pathlib import Path
import numpy as np
from planet_landscape import coupling

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def main():
    path=ROOT/'outputs/review-astra/runs/Chain-360f894f/plates/seed1001/snapshots'
    hist=coupling.PlateHistory(path)
    z=hist.fields(0.)[0]['elevation_m']
    q=hist.g.xyz.copy()
    present=coupling.erodibility_map(hist,q,0.)
    mean=np.zeros(len(q))
    for t in range(0,-30,-1):
        mean+=coupling.erodibility_map(hist,q,float(t))/30.
        q=coupling.step_back(hist,float(t),q)
    rows=[]
    for low,high in [(0,500),(500,2000),(2000,3000),(3000,1e6)]:
        m=(z>low)&(z<=high)
        for tag,v in [('static',present),('parcel_time_mean',mean)]:
            rows.append(dict(set=f'parent_{low}_{high}m',measure=tag,value=np.average(v[m],weights=hist.g.area[m])))
    j=json.loads((ROOT/'outputs/erosion-law-check/868-20261023-203547/results.json').read_text())
    obs=np.array([[r['area_km2'],r['law_m_myr'],r['measured_m_myr']] for r in j['rows']])
    for min_area in (0,500,2000,10000):
        for lo,hi in [(0,20),(20,100),(100,500),(500,np.inf),(0,np.inf)]:
            m=(obs[:,0]>=min_area)&(obs[:,2]>=lo)&(obs[:,2]<hi)
            if m.any():
                key=f'basin_area_ge{min_area}_rate_{lo}_{hi}'
                rows.extend([dict(set=key,measure='n',value=m.sum()),
                             dict(set=key,measure='median_model_over_observed',value=np.median(obs[m,1]/obs[m,2]))])
    with (HERE/'erosion.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    for r in rows: print(r['set'],r['measure'],round(r['value'],3))


if __name__=='__main__':
    main()
