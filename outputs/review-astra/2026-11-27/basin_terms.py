"""One landscape step on observed Earth. Disable each incision mechanism in isolation."""
import csv
import json
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT))
from scripts import erosion_law_check as check
from scripts.run_planet_landscape import THETA
from planet_grid import fields
from planet_grid.sphere_graph import CubeSphere
from planet_landscape import graph_landscape as gl
from planet_earth.bird import lonlat
from matplotlib.path import Path as Polygon
from rasterio.features import rasterize
from rasterio.transform import from_origin


def orogens(g):
    # Same polygon inclusion as planet_earth.bird; cull points by the ring's bounding box.
    lon,lat=lonlat(g.xyz)
    inside=np.zeros(g.size,bool)
    for f in json.loads((ROOT/'data/tectonicplates-master/GeoJSON/PB2002_orogens.json').read_text())['features']:
        geo=f['geometry']
        for poly in geo['coordinates'] if geo['type']=='MultiPolygon' else [geo['coordinates']]:
            ring=np.array(poly[0]);l0=ring[0,0]
            rl=(ring[:,0]-l0+180)%360-180+l0
            pl=(lon-l0+180)%360-180+l0
            m=(pl>=rl.min())&(pl<=rl.max())&(lat>=ring[:,1].min())&(lat<=ring[:,1].max())
            inside[m]|=Polygon(np.c_[rl,ring[:,1]]).contains_points(np.c_[pl[m],lat[m]])
    return inside


def lattice(geom):
    # Rasterize the same exterior-ring, 2 km centre lattice as the archived reader.
    # This avoids a separate point-in-polygon pass over every vertex for every point.
    samples=[]
    for poly in geom['coordinates'] if geom['type']=='MultiPolygon' else [geom['coordinates']]:
        ring=np.array(poly[0])
        (x0,y0),(x1,y1)=ring.min(0),ring.max(0)
        dy=2./111.2
        dx=dy/max(np.cos(np.radians((y0+y1)/2)),.05)
        xs=np.arange(x0,x1+dx,dx);ys=np.arange(y0,y1+dy,dy)
        affine=from_origin(x0-dx/2,ys[-1]+dy/2,dx,dy)
        mask=rasterize([({'type':'Polygon','coordinates':[poly[0]]},1)],out_shape=(len(ys),len(xs)),transform=affine,dtype='uint8')
        iy,ix=np.where(mask)
        samples.append(np.c_[xs[ix],ys[::-1][iy]] if len(ix) else ring[:-1].mean(0,keepdims=True))
    return np.concatenate(samples)


def main():
    g=CubeSphere(868)
    print('Fine geometry ready',flush=True)
    z=np.load(ROOT/'outputs/planet-reference/earth-n868.npz')['elevation_m'].astype(float)
    p={**gl.DEFAULTS,'theta':THETA,'open_edges':False}
    K=p['K']*np.where(orogens(g),1.,.15)
    rates={}
    for tag,k,sc in [('both',K,p['Sc']),('stream',K,np.inf),('slope_cap',K*0,p['Sc'])]:
        out=gl.surface_step(g,z,z< -1500.,None,np.ones(g.size,bool),k,{**p,'Sc':sc},
                            p['dt'],g.spacing_min(),g.diffusion_bound(),keep_deep=True)
        rates[tag]=out['eroded']/p['dt']*1e6
        del out
        print(tag,'ready',flush=True)
    rows=[]
    basins=json.loads((ROOT/'data/octopus/crn_int_basins.json').read_text(encoding='utf-8'))['features']
    for f in basins:
        prop=f['properties']
        obs=check.field(prop,check.RATE_FIELDS)
        area=check.field(prop,check.AREA_FIELDS)
        if not np.isfinite(obs) or obs<=0 or not np.isfinite(area) or area<2000 or not f['geometry']:
            continue
        pts=lattice(f['geometry'])
        ids=fields.nearest_nodes(g,check.to_xyz(pts[:,0],pts[:,1]))
        rows.append(dict(area=area,obs=obs,**{tag:np.mean(v[ids]) for tag,v in rates.items()}))
    with (HERE/'basin_terms_grid.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    for lo,hi in [(0,20),(20,100),(100,500),(500,np.inf),(0,np.inf)]:
        sel=[r for r in rows if lo<=r['obs']<hi]
        print(lo,hi,len(sel),*[round(np.median([r[k]/r['obs'] for r in sel]),3) for k in rates])


if __name__=='__main__':
    main()
