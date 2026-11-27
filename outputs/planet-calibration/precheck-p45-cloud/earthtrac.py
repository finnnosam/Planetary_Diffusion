import sys, numpy as np, earthflow as ef, shx
from planet_reviews import force_balance_earth as fbe, trench_pairs as tp
fn=sys.argv[1]; clip=sys.argv[2]=='slab'; zcut=float(sys.argv[3]) if len(sys.argv)>3 else 300.
lmax,H,V=ef.flow(fn,clip,zcut)
P=ef.spec(H,lmax)
print(f'{fn} {"fast part" if clip else "full"} below {zcut:.0f} km: dyn topo rms {np.sqrt(P[1:].sum()):.0f} m (air-loaded); share l1-3 {P[1:4].sum()/P[1:].sum():.2f}; P30/P2 {P[min(30,lmax)]/P[2]:.4f}; deg rms '+' '.join(f'{l}:{np.sqrt(P[l]):.0f}' for l in (1,2,3,4,6,8,12,20,30)))
g,pid,names,cont=fbe.earth_plates(217); K=len(names)
vfs=shx.grad(V,lmax,g.xyz)                      # m/s
print(f'  free-slip surface flow rms {np.sqrt(np.average((vfs**2).sum(1),weights=g.area))*YR*100:.2f} cm/yr' if (YR:=ef.YR) else '')
st=fbe.earth_steps(names,fbe._rings()); T=fbe.terms(g,pid,cont,st,K); Am=fbe.assemble(T,K); met=fbe.Metric(g,pid,K)
R=fbe.R
F0=np.zeros((K,3)); F1=np.zeros((K,3))
tq=np.cross(g.xyz,vfs)*R*g.area[:,None]
for k in range(K):
    m=pid==k; F0[k]=tq[m].sum(0); F1[k]=(tq[m]*cont[m,None]).sum(0)
Am['F0']=F0.reshape(-1); Am['F1']=F1.reshape(-1)
# coherence of the flow itself on plates (mean flow per plate as rigid fit) skipped
def solve(prm, trac=True):
    Cb,beta,Fsp,Cs,Fsu,Cc=prm
    A=Cb*(Am['D']+beta*Am['K'])+Cs*Am['Z']+Cc*Am['X']
    b=Fsp*Am['S']+fbe.F_RP*Am['P']+Fsu*Am['U']+(Cb*(Am['F0']+beta*Am['F1']) if trac else 0)
    w=np.linalg.solve(A,b); return {k:w[3*k:3*k+3] for k in range(K)}
w_obs=fbe.earth_omega(g,pid,names)
def f(z,s,c):
    try: w=met.nnr(solve(fbe.unpack(z,s,c)))
    except np.linalg.LinAlgError: return 1e9
    return met.misfit(w,w_obs)+fbe.box_penalty(z)
E=tp.earth_segments(); oo=E['io']|E['ocb']; x=st['x']; L=E['L']; cls=E['cls']; ms=R*1e3*ef.YR
def pairs(w,lab):
    va=np.cross(np.stack([w[int(k)] for k in st['A']]),x)*ms; vb=np.cross(np.stack([w[int(k)] for k in st['B']]),x)*ms
    close=-np.einsum('ij,ij->i',vb-va,st['nA']); rho=tp.rho(va,vb)
    row=[f"oo {tp.share(rho[oo]<.5,L[oo]):.2f} (<20: {tp.share(close[oo]<20,L[oo]):.2f})"]+[f"{c} {tp.share(rho[cls==c]<.5,L[cls==c]):.2f}" for c in ('OSR','CRB','OTF','CTF','SUB','OCB','CCB')]
    wn=np.stack([w[int(k)] for k in pid]); coh=tp.coherence(g,pid,np.cross(wn,g.xyz),np.random.default_rng(0))
    print('    rho<0.5: '+'; '.join(row)+'\n    C(d): '+tp.coh_txt(coh))
prm0=fbe.params_coefficients()
for lab,prm in (('PARAMS + traction',prm0),):
    w=met.nnr(solve(prm)); r={'g':g,'pid':pid,'names':names,'st':st,'w':w,'w_obs':w_obs,'K':K,'prm':prm}
    fbe.report(r,'  '+lab); pairs(w,lab)
for s,c in ((False,False),(True,True)):
    z0=np.array([15.,.5,13.5,1.,.2,-1.]); best=None
    for start in (z0,z0+[-.5,-.5,-.3,-.3,.3,1.],z0+[.3,.5,0,.3,0,-1.]):
        z,v=fbe.nelder_mead(lambda q:f(q,s,c),start,[.3,.3,.3,.2,.2,.5]); z,v=fbe.nelder_mead(lambda q:f(q,s,c),z,[.1]*6)
        if best is None or v<best[1]: best=(z,v)
    prm=fbe.unpack(best[0],s,c); w=met.nnr(solve(prm))
    r={'g':g,'pid':pid,'names':names,'st':st,'w':w,'w_obs':w_obs,'K':K,'prm':prm}
    fbe.report(r,f'  refit with traction, suction {s}, coupling {c}  (box z={np.round(best[0],2)})'); pairs(w,'')
    # flow off, same coefficients
    w2=met.nnr(solve(prm,False)); print('   same coefficients, traction off:'); pairs(w2,'')
