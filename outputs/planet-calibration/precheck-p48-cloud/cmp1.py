import numpy as np, dmf, colloc
a=6371e3; c=dmf.CMB_RADIUS_M
# 1) isoviscous: both solvers
base=np.array([c/a]); eta=np.array([1e21])
zs=np.array([325,600,1000,1500,2000,2500,2825])*1e3
for N in (16,24,32):
    pass
K=dmf.FlowKernels(a,a-zs,base,eta,31)
print('isoviscous, KH (D | colloc N=24 | N=40)')
for l in (1,2,3,8,20,31):
    row=[]
    for j,z in enumerate(zs):
        kh,kv=colloc.kernel(l,(a-z)/a,[c/a,1.],[1e21],N=24)
        kh2,kv2=colloc.kernel(l,(a-z)/a,[c/a,1.],[1e21],N=40)
        row.append(f"{K.KH[l,j]:+.5f}|{kh:+.5f}|{kh2:+.5f}")
    print(l,'  '.join(row))
print('isoviscous KV ratio colloc/D')
for l in (1,2,8,31):
    print(l,[round(colloc.kernel(l,(a-z)/a,[c/a,1.],[1e21],N=32)[1]/K.KV[l,j],6) for j,z in enumerate(zs)])
