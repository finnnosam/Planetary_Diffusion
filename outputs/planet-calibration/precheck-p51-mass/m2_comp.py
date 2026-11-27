"""Pass 51, M2: compositional mass of subducted crust + depleted residue in the lower mantle (P51_RULE.txt), per m^2
of plate, from normative lower-mantle modes and SLB2011 volumes on the ambient isentrope. argv: [K_D] [crust km]
[column km]. Writes m2_<KD>_<crust>_<col>.json"""
import json, sys
import numpy as np
import slb_common as s
from burnman.minerals import SLB_2011 as S

KD = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
HC = float(sys.argv[2]) if len(sys.argv) > 2 else 7.1
HCOL = float(sys.argv[3]) if len(sys.argv) > 3 else 60.
MW = dict(SiO2=60.08, Al2O3=101.96, FeO=71.85, MgO=40.30, CaO=56.08, Na2O=61.98)
PYR = dict(SiO2=45.0, Al2O3=4.45, FeO=8.05, MgO=37.8, CaO=3.55, Na2O=0.36)        # McDonough & Sun 1995
MORB = dict(SiO2=50.06, Al2O3=14.70, FeO=10.43, MgO=7.58, CaO=11.39, Na2O=2.79)   # Gale et al. 2013 all-MORB mean
EM = dict(mgpv=S.mg_perovskite(), fepv=S.fe_perovskite(), alpv=S.al_perovskite(), pe=S.periclase(), wu=S.wuestite(),
          capv=S.ca_perovskite(), st=S.stishovite(), nacf=S.na_ca_ferrite())


def norm(wt):
    tot = sum(wt.values()); return {k: v/tot*100 for k, v in wt.items()}


def modes(wt):
    """mol of each endmember per 100 g (normative lower-mantle allocation)."""
    m = {k: wt[k]/MW[k] for k in MW}
    Si, Al, Fe, Mg, Ca, Na = m['SiO2'], 2*m['Al2O3'], m['FeO'], m['MgO'], m['CaO'], 2*m['Na2O']
    out = dict(capv=Ca, nacf=Na); Si -= Ca+Na; Al -= Na
    out['alpv'] = Al/2
    nb = min(Si, Mg+Fe); out['st'] = max(Si-nb, 0.); nf = Mg+Fe-nb
    if nf <= 1e-12:
        out.update(mgpv=Mg, fepv=Fe, pe=0., wu=0.)
    else:   # split Fe between bdg (nb) and fp (nf) with (xb/(1-xb))/(xf/(1-xf)) = KD
        lo, hi = 0., 1.
        for _ in range(80):
            xf = 0.5*(lo+hi); rf = xf/(1-xf); xb = KD*rf/(1+KD*rf)
            if xb*nb+xf*nf > Fe: hi = xf
            else: lo = xf
        out.update(mgpv=(1-xb)*nb, fepv=xb*nb, pe=(1-xf)*nf, wu=xf*nf)
    return out


def specific_volume(wt, P, T):
    v = 0.; mass = 0.
    for k, n in modes(wt).items():
        e = EM[k]; e.set_state(P, T)
        v += n*e.molar_volume; mass += n*e.molar_mass
    return v/mass                                            # m^3/kg


pyr, morb = norm(PYR), norm(MORB)
m_c = HC*1e3*2900.                                           # crust mass per m^2 of plate (surface density 2,900)
m_col = HCOL*1e3*3300.
res = {k: (pyr[k]*m_col-morb[k]*m_c)/(m_col-m_c) for k in pyr}   # residue (wt%): the column minus its crust
ZL = 300+25+50*np.arange(51); Ta = s.isentrope(1600., 0., ZL)
rows = []
for z, T in zip(ZL, Ta):
    if z < 670:
        rows.append(0.); continue
    P = s.P_of_z(z); vp = specific_volume(pyr, P, T)
    vc = specific_volume(morb, P, T); vr = specific_volume(res, P, T)
    dm = m_c*(1-vc/vp)+(m_col-m_c)*(1-vr/vp)                # sum over the column of (1 - rho_pyr/rho) dm
    rows.append(dm)
    if z in (725, 1025, 1525, 2025, 2525, 2825):
        print(f'{z} km: MORB {100*(vp/vc-1):+.2f}%  residue {100*(vp/vr-1):+.3f}%  crust part {m_c*(1-vc/vp):+.3e}  residue part '
              f'{(m_col-m_c)*(1-vr/vp):+.3e}  net {dm:+.3e} kg/m2')
md = modes(pyr); tot = sum(n*EM[k].molar_mass for k, n in md.items())
print('pyrolite modes (mass):', {k: round(n*EM[k].molar_mass/tot, 3) for k, n in md.items()})
mm = modes(morb); tot = sum(n*EM[k].molar_mass for k, n in mm.items())
print('MORB modes (mass):', {k: round(n*EM[k].molar_mass/tot, 3) for k, n in mm.items()})
json.dump(dict(KD=KD, crust_km=HC, column_km=HCOL, z_km=ZL.tolist(), dm_kg_m2=rows, residue_wt=res),
          open(f'm2_{KD}_{HC}_{HCOL:.0f}.json', 'w'), indent=1)
