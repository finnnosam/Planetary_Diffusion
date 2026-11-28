"""Run from any directory with the project's Python. Writes one row per map."""
import csv
import json
from pathlib import Path

import numpy as np
from planet_earth.reference import ice_mask
from planet_earth.crust1 import OCEANIC, crust1_cell
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics.graphops import components, nearest_label

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CAL = ROOT / "outputs/planet-calibration"
REF = ROOT / "outputs/planet-reference"
G = CubeSphere(217)
A = G.area
I, K = np.where(G.axis & (G.nbr > np.arange(G.size)[:, None]))
J = G.nbr[I, K]
DIST = G.dist[I, K]
EDGE_W = np.sqrt(A[I] * A[J]) / DIST


def quant(v, w, q=.5):
    order = np.argsort(v)
    cum = np.cumsum(w[order])
    return v[order][np.searchsorted(cum, np.asarray(q) * cum[-1])]


def assess(name, z, valid, thickness, continental, age, age_w, plate):
    land = (z > 0) & valid
    wet = (z <= 0) & valid
    def share(mask, domain):
        return 100 * A[mask & domain].sum() / A[domain].sum()
    r = dict(name=name, land_pct=share(land, valid), land_mean_m=np.average(z[land], weights=A[land]))
    for p, v in zip([10, 50, 90, 99], quant(z[land], A[land], [.1,.5,.9,.99])):
        r[f'land_p{p}']=v
    for p, v in zip([10, 50, 90], quant(z[wet], A[wet], [.1,.5,.9])):
        r[f'ocean_p{p}']=v
    for tag, mask in [('low',z<500), ('mid',(z>=500)&(z<2000)),('high3',z>3000),('high4',z>4000)]:
        r[tag+'_land_pct']=share(mask, land)
    for tag, mask in [('land',land),('top1',land & (z>=r['land_p99'])),('high3',land&(z>3000)),('low',land&(z<500))]:
        labels = components(mask, G.nbr)
        size = np.sort(np.bincount(labels[mask], weights=A[mask]))[::-1]/1e12
        r[tag+'_largest_Mkm2'] = size[0] if size.size else 0
        r[tag+'_largest_pct'] = 100*size[0]/size.sum() if size.size else 0
    # Distance to an actual wet cell. Ignored polar cells can neither seed nor create an artificial coast.
    _, coast_dist = nearest_label(np.where(z<=0, 0, -1), G.nbr, G.dist)
    r['inland_500_pct'] = share(coast_dist>500e3,land)
    r['inland_1000_pct'] = share(coast_dist>1e6,land)
    r['wet_0_200_pct'] = share(z>=-200,wet)
    r['wet_200_2000_pct']=share((z< -200)&(z>=-2000),wet)
    r['thick55_land_pct']=share(continental & (thickness>55000),land)
    for lo,hi in [(30,35),(35,40),(40,45),(45,55),(55,80)]:
        m=valid & continental & (thickness>=lo*1000)&(thickness<hi*1000)
        r[f'z_T{lo}_{hi}']=quant(z[m],A[m]) if m.any() else np.nan
    for tag, mask in [('land',land),('deep',wet&(z< -3000)),('inland',land&(coast_dist>500e3)),('low',land&(z<500))]:
        edge=mask[I]&mask[J]
        r[tag+'_edge_med_m']=quant(abs(z[I[edge]]-z[J[edge]]),EDGE_W[edge])
    coast=(land[I]&wet[J])|(wet[I]&land[J])
    r['coast_proxy_km']=EDGE_W[coast].sum()/1e3
    r['coast_step_m']=quant(abs(z[I[coast]]-z[J[coast]]),EDGE_W[coast])
    # Shared cell-face crossings, measured identically for the plate map and PB2002.
    boundary=(plate[I]!=plate[J])&(plate[I]>=0)&(plate[J]>=0)
    r['plate_edge_proxy_km']=EDGE_W[boundary].sum()/1e3
    good=np.isfinite(age)&(age_w>0)
    r['age_mean']=np.average(age[good],weights=(age_w*A)[good])
    r['age_lt10_Mkm2']=np.sum((age_w*A)[good&(age<10)])/1e12
    return r


def main():
    z=np.load(REF/'earth-n217.npz')['elevation_m']
    t=np.loadtxt(ROOT/'data/crust1/crsthk.xyz')[:,2][crust1_cell(G)]*1000
    types=np.loadtxt(ROOT/'data/crust1/CNtype1-1.txt',dtype=str).ravel()[crust1_cell(G)]
    ag=np.load(REF/'earth-seafloor-age-n217.npz')
    plates=np.load(REF/'earth-plates-n217.npz')['plate']
    rows=[assess(name,z,mask,t,~np.isin(types,list(OCEANIC)),ag['age_myr'],ag['dated_share'],plates)
          for name,mask in [('Earth',~ice_mask(G)),('EarthAll',np.ones(G.size,bool))]]
    sources=[('nearest','sweep-20261126-161947/remap_rule-nearest_rift_area-world'),
             ('overlap','sweep-20261126-141105/remap_rule-overlap'),
             ('overlap_plates','sweep-20261126-161947/remap_rule-overlap_rift_area-plates')]
    for arm,folder in sources:
        for seed in sorted((CAL/folder).glob('seed*')):
            s=json.loads((seed/'stats.json').read_text())
            with np.load(seed/'snapshots/t+0000.0.npz') as f:
                z=f['elevation_m']
                c=f['continental'].astype(bool)
                r=assess(f'{arm}/{seed.name}',z,np.ones(G.size,bool),f['crust_thickness_m']+f['sediment_bulk_m'],
                         c,f['ocean_age_myr'],~c,f['plate_id'])
            assert abs(r['land_pct']/100-s['present']['land_fraction']) < 1e-5
            rows.append(r)
    with (HERE/'spatial.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    for key in list(rows[0])[1:]:
        print(f'{key:26}',*[f'{row[key]:9.2f}' for row in rows[:2]], *[f'{np.median([r[key] for r in rows if r["name"].startswith(arm+"/")]):9.2f}' for arm,_ in sources])


if __name__=='__main__':
    main()
