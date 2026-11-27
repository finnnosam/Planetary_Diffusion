"""Pass 51, M1: the slab's thermal mass at depth from SLB2011 (P51_RULE.txt), per unit of its surface (GDH1) mass,
against D's alpha(z)/alpha_GDH1. Plate: GDH1 (Stein & Stein 1992: a 95 km, T_m 1450 C, kappa 0.804e-6) deficit
T_m - T(y) laid on the ambient isentrope; each element follows its own isentrope; no-conduction and fully diffused
limits. argv: [T_pot K] [age Myr]. Writes m1_<Tpot>_<age>.json"""
import json, sys
import numpy as np
import slb_common as s
from planet_tectonics import mantle_flow as mf

TP = float(sys.argv[1]) if len(sys.argv) > 1 else 1600.
AGE = float(sys.argv[2]) if len(sys.argv) > 2 else 63.
A_PL, TM, KAP = 95e3, 1450., 0.804e-6


def gdh1_deficit(y, age):
    t = age*3.15576e13; th = y/A_PL
    T = th.copy()
    for n in range(1, 200):
        T += 2/(n*np.pi)*np.sin(n*np.pi*th)*np.exp(-n*n*np.pi**2*KAP*t/A_PL**2)
    return TM*(1-T)


ZL = 300+25+50*np.arange(51)                       # D's load-layer mids (km), 325 ... 2825
y = np.arange(1., 95., 2.)*1e3                     # element mid-depths in the plate (m), 2 km thick
dT = gdh1_deficit(y, AGE)
Ta_y = s.isentrope(TP, 0., y/1e3)
Te_y = Ta_y-dT
# reference: the plate column at the trench (upper-mantle rock at the element's own depth)
R0 = s.rock(0)
ref = 0.; dm = []
for yy, ta, te in zip(y, Ta_y, Te_y):
    P = s.P_of_z(yy/1e3)
    ra = R0.props(P, ta)[0]; re = R0.props(P, te)[0]
    dm.append(re*2e3); ref += (1-ra/re)*re*2e3
dm = np.array(dm)
Ta_z = s.isentrope(TP, 0., ZL)
Te = np.array([s.isentrope(te, yy/1e3, ZL) for yy, te in zip(y, Te_y)])     # (elements, layers)
F_nc, F_df, F_lin = [], [], []
for j, z in enumerate(ZL):
    reg = s.region(z); R = s.rock(reg); P = s.P_of_z(z)
    ra, aa, cpa = R.props(P, Ta_z[j])
    m_nc = 0.; Q = 0.
    for i in range(len(y)):
        re, ae, cpe = R.props(P, Te[i, j])
        m_nc += (1-ra/re)*dm[i]; Q += 0.5*(cpa+cpe)*(Ta_z[j]-Te[i, j])*dm[i]
    F_nc.append(m_nc/ref); F_df.append(aa*Q/cpa/ref)
    F_lin.append(aa/mf.ALPHA_GDH1*(Ta_z[j]-Te[:, j]).dot(dm)/dT.dot(dm))
F_D = mf.alpha_ratio(ZL)
out = dict(T_pot=TP, age=AGE, z_km=ZL.tolist(), F_noconduction=F_nc, F_diffused=F_df, F_D=F_D.tolist(),
           ref_kg_m2=ref, gdh1_kg_m2=float(mf.thermal_mass_excess(np.array([AGE]), 3300., 1030.)[0]),
           deficit_amplification=((Ta_z[None, :]-Te).T.dot(dm)/dT.dot(dm)).tolist(), Ta_z=Ta_z.tolist())
json.dump(out, open(f'm1_{TP:.0f}_{AGE:.0f}.json', 'w'), indent=1)
print(f'T_pot {TP:.0f} K, age {AGE:.0f} Myr: plate mass SLB {ref:.3e} vs GDH1 {out["gdh1_kg_m2"]:.3e} kg/m2')
for j in range(0, len(ZL), 5):
    print(f'{ZL[j]:5.0f} km: F no-cond {F_nc[j]:.3f}  diffused {F_df[j]:.3f}  D {F_D[j]:.3f}  ratio {0.5*(F_nc[j]+F_df[j])/F_D[j]:.2f}  '
          f'deficit x{out["deficit_amplification"][j]:.2f}')
