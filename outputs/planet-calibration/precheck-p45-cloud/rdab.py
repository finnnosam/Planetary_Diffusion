import numpy as np, ducc0
def read_ab(fn):
    tok=open(fn).read().split()
    nl=int(tok[0]); i=1; out=[]
    for _ in range(nl):
        z=float(tok[i]); lmax=int(tok[i+1]); i+=2; n=(lmax+1)*(lmax+2)//2
        v=np.array(tok[i:i+2*n],float).reshape(n,2); i+=2*n; out.append((z,v))
    out.sort(key=lambda t:t[0]); return lmax,out
def to_alm(v,lmax):
    a=np.zeros((lmax+1)*(lmax+2)//2,complex); k=0
    for l in range(lmax+1):
        for m in range(l+1):
            A,B=v[k]; k+=1; a[m*(2*lmax+1-m)//2+l]=A if m==0 else (A-1j*B)/np.sqrt(2)
    return a
