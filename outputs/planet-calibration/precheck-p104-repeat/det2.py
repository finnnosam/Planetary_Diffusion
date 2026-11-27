import numpy as np, ducc0, time
rng = np.random.default_rng(1)
def loc_of(xyz):
    return np.c_[np.arccos(np.clip(xyz[:,2],-1,1)), np.mod(np.arctan2(xyz[:,1], xyz[:,0]), 2*np.pi)]
# clustered: points along a few great-circle girdles with small scatter (slab parcels)
pts = []
for g in range(6):
    a = rng.normal(size=3); a/=np.linalg.norm(a); b = np.cross(a, rng.normal(size=3)); b/=np.linalg.norm(b)
    t = rng.uniform(0, 2*np.pi, 30000)
    p = np.outer(np.cos(t), a)+np.outer(np.sin(t), b)+rng.normal(scale=.01, size=(30000,3))
    pts.append(p/np.linalg.norm(p, axis=1)[:,None])
xyz = np.vstack(pts); loc = loc_of(xyz); w = rng.uniform(size=len(xyz))[None]
for lmax, eps in ((30, 1e-10), (480, 1e-7)):
    for nt in (1, 2):
        ref = ducc0.sht.experimental.adjoint_synthesis_general(map=w, spin=0, lmax=lmax, loc=loc, epsilon=eps, nthreads=nt)
        nd = 0; t0 = time.time()
        for k in range(60):
            o = ducc0.sht.experimental.adjoint_synthesis_general(map=w, spin=0, lmax=lmax, loc=loc, epsilon=eps, nthreads=nt)
            nd += not np.array_equal(o, ref)
        print(f'lmax {lmax} nthreads {nt}: {nd} of 60 repeats differ; {(time.time()-t0)/60*1e3:.1f} ms each')
