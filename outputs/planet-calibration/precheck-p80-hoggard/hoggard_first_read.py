# First read of Hoggard et al. 2020 SL2013sv lithospheric thickness on CRUST1.0 continental cells (9 Nov, after pass 79):
# L by tectonothermal age class (lithosphere.CLASS_AGE_MA) and by crust thickness. Run from the repo root.
import numpy as np, sys
sys.path.insert(0,'packages/planet-tectonics/src')
from planet_tectonics.lithosphere import CLASS_AGE_MA
H='data/Hoggard2020/twksd-osfstorage-archive/supplementary_datasets/lithospheric_thickness_maps/global/SL2013sv/SL2013sv.csv'
d=np.loadtxt(H,delimiter=','); lo,la,v=d.T
grid=np.full((361,721),np.nan); grid[np.round((90-la)*2).astype(int),np.round(lo*2).astype(int)]=v
types=np.array(open('data/crust1/CNtype1-1.txt').read().split()).reshape(180,360)
th=np.loadtxt('data/crust1/crsthk.xyz')[:,2].reshape(180,360)
lat=np.arange(89.5,-90,-1); lon=np.arange(-179.5,180,1)
LA,LO=np.meshgrid(lat,lon,indexing='ij')
L=grid[np.round((90-LA)*2).astype(int),np.round((LO%360)*2).astype(int)]
w=np.cos(np.radians(LA))
age=np.vectorize(lambda t: CLASS_AGE_MA.get(t,np.nan))(types)
c=np.isfinite(age)&(th>=20)
print('continental CRUST1.0 cells with a class age:',c.sum())
print('age class (Ma) | area share | Hoggard SL2013sv L p25/p50/p75 | C law 93.6+0.04t')
for a in sorted(set(age[c])):
    m=c&(age==a); print(f'{a:6.0f} | {w[m].sum()/w[c].sum():.3f} | {np.round(np.percentile(L[m],[25,50,75]))} | {93.6+.04*a:.0f}')
x=age[c]; y=L[c]; ww=w[c]
A=np.vstack([np.ones_like(x),x]).T*np.sqrt(ww)[:,None]; b=np.linalg.lstsq(A,y*np.sqrt(ww),rcond=None)[0]
print('weighted least squares L = %.1f + %.4f t (km, Ma); C uses 93.6 + 0.0400 t'%tuple(b))
print('crust T bin | L p50 (all continental types)')
for a_,b_ in [(20,25),(25,30),(30,35),(35,40),(40,45),(45,55),(55,80)]:
    m=c&(th>=a_)&(th<b_); print(f'  {a_}-{b_} | {np.median(L[m]):.0f} (area {w[m].sum()/w[c].sum():.3f})')
print('dL/dT (weighted) %.2f km per km'%np.polyfit(th[c],L[c],1,w=np.sqrt(w[c]))[0])
