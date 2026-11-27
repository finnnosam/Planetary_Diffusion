import numpy as np
from planet_tectonics import mantle_flow as mf
from planet_grid.paths import DATA
A=6371e3; E=np.load('e3_loads.npz'); RAD=E['radii']; L=int(E['lmax']); FAC=mf.alpha_ratio((A-RAD)/1e3)
base,eta=mf.read_viscosity(DATA/'mantle'/'visc.sh08'); K=mf.FlowKernels(A,RAD,base,eta,L)
loads=E['loads_0']*FAC[:,None]
H,_=K.surface(loads); tot=np.sqrt(mf.spectrum(H/3300.,L)[1:].sum()); print('total',tot)
z=(A-RAD)/1e3
bins=[(300,660),(660,1000),(1000,1500),(1500,2000),(2000,2500),(2500,2850)]
mb=E['mass_by_layer_0']*FAC
for a,b in bins:
    s=(z>a)&(z<b); l2=loads.copy(); l2[~s]=0
    Hb,_=K.surface(l2); r=np.sqrt(mf.spectrum(Hb/3300.,L)[1:].sum())
    # projection onto total
    li=mf.alm_l(L); W=np.where(mf.alm_m(L)==0,1.,2.); ok=li>=1
    proj=(W[ok]*(Hb[ok]*np.conj(H[ok])).real).sum()/(W[ok]*abs(H[ok])**2).sum()
    print(f'{a}-{b}: alone {r:.0f} m, share of total power (projection) {proj:.2f}, mass*alpha {mb[s].sum()/1e21:.2f}e21')
# kernel KH at l=2 by depth
for zz in (350,500,700,1000,1500,2000,2500,2800):
    j=np.argmin(abs(z-zz)); print(zz, 'KH l1..4', np.round(K.KH[1:5,j],3), 'KN?')
