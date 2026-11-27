import numpy as np
A_M=6371e3; C_M=3480e3; G=9.81
def amat(l,e):
    L=l*(l+1.)
    return np.array([[-2,L,0,0],[-1,1,0,1/e],[12*e,-6*L*e,1,L],[-6*e,2*(2*L-1)*e,-1,-2]],float)
def prop_eig(l,e,x):
    w,W=np.linalg.eig(amat(l,e)); w=w.real; W=W.real
    return (W*np.exp(w*x))@np.linalg.inv(W)
def prop_rk(l,e,x,n=4000):
    A=amat(l,e); P=np.eye(4); h=x/n
    for _ in range(n):
        k1=A@P; k2=A@(P+h/2*k1); k3=A@(P+h/2*k2); k4=A@(P+h*k3); P=P+h/6*(k1+2*k2+2*k3+k4)
    return P
def P_between(l,r2,r1,rb,eta,eta0,f=prop_eig):
    # rb: layer top radii ascending (last = A_M), eta: layer viscosity; propagate r1->r2 (r2>r1)
    P=np.eye(4); r=r1
    for top,et in zip(rb,eta):
        if top<=r: continue
        nxt=min(top,r2)
        P=f(l,et/eta0,np.log(nxt/r))@P; r=nxt
        if r>=r2: break
    return P
def solve(l,ri,rb,eta,eta0=1e21,f=prop_eig):
    """unit surface density (1 kg/m2) at radius ri. returns U,V(a), sigma_rr(a), sigma_rr(c)"""
    Pac=P_between(l,A_M,C_M,rb,eta,eta0,f); Pai=P_between(l,A_M,ri,rb,eta,eta0,f)
    b=np.array([0,0,ri*G*1.0/eta0,0])
    s=Pai@b
    # X(a)=Pac@[0,V,Rc,0]+s ; need X0(a)=0, X3(a)=0
    M=Pac[np.ix_([0,3],[1,2])]; rhs=-s[[0,3]]
    Vc,Rc=np.linalg.solve(M,rhs)
    Xa=Pac@np.array([0,Vc,Rc,0])+s
    return Xa, np.array([0,Vc,Rc,0])
