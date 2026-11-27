import numpy as np, ducc0
def idx(l,m,lmax): return m*(2*lmax+1-m)//2+l
def Lsqrt(lmax):
    out=np.zeros((lmax+1)*(lmax+2)//2)
    for m in range(lmax+1):
        for l in range(m,lmax+1): out[idx(l,m,lmax)]=np.sqrt(l*(l+1.))
    return out
def loc_of(xyz):
    th=np.arccos(np.clip(xyz[:,2],-1,1)); ph=np.arctan2(xyz[:,1],xyz[:,0])%(2*np.pi)
    return np.stack([th,ph],1)
def syn(a,lmax,xyz):
    return ducc0.sht.experimental.synthesis_general(alm=a[None],spin=0,lmax=lmax,loc=loc_of(xyz),epsilon=1e-10)[0]
def grad(a,lmax,xyz,sign=1.):
    """surface gradient (unit sphere) of the scalar with coefficients a, as xyz vectors"""
    loc=loc_of(xyz)
    g=ducc0.sht.experimental.synthesis_general(alm=np.stack([sign*a*Lsqrt(lmax),0*a]),spin=1,lmax=lmax,loc=loc,epsilon=1e-10)
    th,ph=loc[:,0],loc[:,1]
    eth=np.stack([np.cos(th)*np.cos(ph),np.cos(th)*np.sin(ph),-np.sin(th)],1)
    eph=np.stack([-np.sin(ph),np.cos(ph),0*ph],1)
    return g[0][:,None]*eth+g[1][:,None]*eph
