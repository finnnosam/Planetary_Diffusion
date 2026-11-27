"""Pass 45 pre-check: slab-driven free-slip flow from Earth's tomography (fast part below Z_CUT), its dynamic topography
spectrum, and traction toward it in the torque law on Bird's plates."""
import sys, numpy as np, mf, shx, rdab
YR=3.15576e7
PZ=np.array([0,220,400,670,771,1000,1500,2000,2500,2891.]); PR=np.array([3380,3360,3540,3990,4440,4580,4800,5050,5300,5560.])
def visc(fn='/home/claude/d/visc.sh08'):
    v=np.loadtxt(fn); return list(v[1:,0]*mf.A_M)+[mf.A_M], list(v[:,1])
def kernels(lmax,depths_km,rb,eta):
    KH=np.zeros((lmax+1,len(depths_km))); KV=np.zeros_like(KH)
    for j,z in enumerate(depths_km):
        ri=mf.A_M-z*1e3
        for l in range(1,lmax+1):
            Xa,_=mf.solve(l,ri,rb,eta)
            KH[l,j]=-(Xa[2]*1e21/mf.A_M)/mf.G      # surface normal stress / g  (kg/m2) per kg/m2 load ; h = KH*sig/drho
            KV[l,j]=Xa[1]                           # m/s per kg/m2
    return KH,KV
def density_loads(fn,clip,zcut,dfac=.2):
    lmax,L=rdab.read_ab(fn); zs=np.array([z for z,_ in L]); edges=np.concatenate([[0],(zs[1:]+zs[:-1])/2,[2891]])
    nth,nph=2*lmax+2,4*lmax+4; out=[]
    import ducc0
    for k,(z,v) in enumerate(L):
        if z<zcut: continue
        a=rdab.to_alm(v,lmax)*dfac*0.01*np.interp(z,PZ,PR)
        if clip:
            m=ducc0.sht.synthesis_2d(alm=a[None],spin=0,lmax=lmax,geometry='GL',ntheta=nth,nphi=nph)[0]
            a=ducc0.sht.analysis_2d(map=np.maximum(m,0)[None],spin=0,lmax=lmax,geometry='GL')[0]
        dz=(edges[k+1]-max(edges[k],zcut))*1e3
        out.append((z,a*dz))
    return lmax,out
def flow(fn,clip,zcut=300.,drho=3300.):
    rb,eta=visc(); lmax,loads=density_loads(fn,clip,zcut)
    KH,KV=kernels(lmax,[z for z,_ in loads],rb,eta)
    n=(lmax+1)*(lmax+2)//2; H=np.zeros(n,complex); V=np.zeros(n,complex)
    lof=np.array([l for m in range(lmax+1) for l in range(m,lmax+1)])
    for j,(z,s) in enumerate(loads):
        H+=-KH[lof,j]*s/drho; V+=KV[lof,j]*s
    return lmax,H,V
def spec(a,lmax):
    P=np.zeros(lmax+1)
    for m in range(lmax+1):
        for l in range(m,lmax+1): P[l]+=abs(a[shx.idx(l,m,lmax)])**2*(1 if m==0 else 2)
    return P/(4*np.pi)
