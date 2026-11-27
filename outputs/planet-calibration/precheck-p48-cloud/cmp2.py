import numpy as np, dmf, colloc
a=6371e3; c=dmf.CMB_RADIUS_M
base,eta=dmf.read_viscosity('visc.sh08')
edges=np.append(np.append(c/a, base[1:]),1.)
radii=dmf.load_radii(a,300e3)
K=dmf.FlowKernels(a,radii,base,eta,31)
worst_h=0; worst_v=0
out=[]
for l in range(1,32):
    for j,r in enumerate(radii):
        kh,kv=colloc.kernel(l,r/a,edges,eta,N=16)
        dh=abs(kh-K.KH[l,j]); worst_h=max(worst_h,dh)
        if abs(K.KV[l,j])>1e-3*abs(K.KV[l]).max(): worst_v=max(worst_v,abs(kv/K.KV[l,j]-1))
        out.append((l,(a-r)/1e3,K.KH[l,j],kh,K.KV[l,j],kv))
print('max |KH_D - KH_ref| over l 1-31, 51 load depths:',worst_h)
print('max rel diff KV (where |KV| > 1e-3 of its degree max):',worst_v)
for l,z,a1,b1,c1,d1 in out:
    if l in (1,2,3,8,20) and int(z) in (325,575,1025,1525,2025,2525,2825):
        print(f"l={l:2d} z={z:5.0f}  KH D {a1:+.5f} ref {b1:+.5f}   KV D {c1:+.4e} ref {d1:+.4e}")
np.save('kernels_ref.npy',np.array(out))
