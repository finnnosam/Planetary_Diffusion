import numpy as np, dmf, colloc
a=6371e3; c=dmf.CMB_RADIUS_M
base,eta=dmf.read_viscosity('visc.sh08')
edges=np.append(np.append(c/a, base[1:]),1.)
radii=dmf.load_radii(a,300e3)
K=dmf.FlowKernels(a,radii,base,eta,31)
for l,z in ((1,1925),(1,1775),(2,1925),(1,1375),(1,525)):
    j=int(np.argmin(abs((a-radii)/1e3-z)))
    vals=[colloc.kernel(l,radii[j]/a,edges,eta,N=N)[0] for N in (12,16,20,28,36)]
    print(l,z,'D',round(K.KH[l,j],5),'colloc N=12..36',np.round(vals,5))
