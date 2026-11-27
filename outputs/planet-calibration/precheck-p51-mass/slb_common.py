"""Pass 51 shared pieces: PREM pressure, SLB2011 pyrolite composites by region, isentropes (dT/dP = alpha T/(rho c_p),
measured latent-heat jumps at 410 and 660 applied as fractional jumps)."""
import warnings
import numpy as np
warnings.filterwarnings('ignore')
import burnman
from burnman.minerals import SLB_2011 as S

G = 6.674e-11; RE = 6371e3


def prem_rho(r):
    x = r/RE; z = RE-r
    if r < 1221.5e3: return 1e3*(13.0885-8.8381*x**2)
    if r < 3480e3: return 1e3*(12.5815-1.2638*x-3.6426*x**2-5.5281*x**3)
    if z > 670e3: return 1e3*(7.9565-6.4761*x+5.5283*x**2-3.0807*x**3)
    if z > 600e3: return 1e3*(5.3197-1.4836*x)
    if z > 400e3: return 1e3*(11.2494-8.0298*x)
    if z > 220e3: return 1e3*(7.1089-3.8045*x)
    if z > 24.4e3: return 1e3*(2.6910+0.6924*x)
    return 2900.


_r = np.linspace(0, RE, 64001); _rho = np.array([prem_rho(r) for r in _r])
_m = np.concatenate([[0], np.cumsum(4*np.pi*((_r[1:]+_r[:-1])/2)**2*(_rho[1:]+_rho[:-1])/2*np.diff(_r))])
_g = np.where(_r > 0, G*_m/np.maximum(_r, 1)**2, 0)
_dP = (_rho*_g)[::-1]; _rr = _r[::-1]
_P = np.concatenate([[0], np.cumsum((_dP[1:]+_dP[:-1])/2*(-np.diff(_rr)))])
ZTAB = (RE-_rr)/1e3; PTAB = _P


def P_of_z(z_km):
    return np.interp(z_km, ZTAB, PTAB)


def _ss(a, b, f):  # binary solid solution as endmembers in a composite
    return [(a, 1-f), (b, f)]


def region(z_km):
    return 0 if z_km < 410 else 1 if z_km < 520 else 2 if z_km < 670 else 3


def composite(reg, fe=0.10):
    """Pyrolite modes (approximate volume fractions; SLB2011 endmembers). 0: ol+opx+cpx+gt; 1: wad+gt; 2: rw+gt+capv;
    3: bdg+fp+capv."""
    if reg == 0:
        parts = [(S.forsterite(), .60*(1-fe)), (S.fayalite(), .60*fe), (S.enstatite(), .15), (S.diopside(), .12),
                 (S.pyrope(), .13*.8), (S.almandine(), .13*.2)]
    elif reg == 1:
        parts = [(S.mg_wadsleyite(), .60*(1-fe)), (S.fe_wadsleyite(), .60*fe), (S.pyrope(), .40*.7),
                 (S.mg_majorite(), .40*.2), (S.grossular(), .40*.1)]
    elif reg == 2:
        parts = [(S.mg_ringwoodite(), .60*(1-fe)), (S.fe_ringwoodite(), .60*fe), (S.pyrope(), .33*.75),
                 (S.mg_majorite(), .33*.25), (S.ca_perovskite(), .07)]
    else:
        parts = [(S.mg_perovskite(), .75*.88), (S.fe_perovskite(), .75*.05), (S.al_perovskite(), .75*.07),
                 (S.periclase(), .18*.85), (S.wuestite(), .18*.15), (S.ca_perovskite(), .07)]
    return parts


class Rock:
    """Volume-fraction mixture: rho = sum f_i rho_i at the reference state's fractions held as mass fractions."""
    def __init__(self, parts, Pref=None, Tref=1600.):
        Pref = Pref if Pref is not None else 1e5
        self.min = [m for m, _ in parts]
        f = np.array([v for _, v in parts]); f /= f.sum()
        rho = []
        for m in self.min:
            m.set_state(Pref, Tref); rho.append(m.density)
        w = f*np.array(rho); self.w = w/w.sum()            # mass fractions (fixed)

    def props(self, P, T):
        v = 0.; a = 0.; cp = 0.
        for m, w in zip(self.min, self.w):
            m.set_state(P, T)
            vi = w/m.density; v += vi; a += vi*m.alpha
            cp += w*m.molar_heat_capacity_p/m.molar_mass
        return 1/v, a/v, cp                                 # density, alpha (volume-weighted), c_p per kg


ROCKS = {}


def rock(reg):
    if reg not in ROCKS:
        Pref = {0: 1e5, 1: 14e9, 2: 19e9, 3: 25e9}[reg]
        ROCKS[reg] = Rock(composite(reg), Pref)
    return ROCKS[reg]


JUMP = {410: +43./1750., 660: -20./2000.}   # fractional T jumps of the isentrope (Clapeyron 3 / -1 MPa/K; see the rule)


def isentrope(T0, z0_km, z_km, dz=5.):
    """T along the isentrope from (z0, T0) to each depth in z_km (increasing)."""
    out = []; z = z0_km; T = T0
    for zt in z_km:
        while z < zt-1e-9:
            h = min(dz, zt-z); reg = region(z+h/2)
            P0, P1 = P_of_z(z), P_of_z(z+h)
            rho, a, cp = rock(reg).props((P0+P1)/2, T)
            T = T*np.exp(a/(rho*cp)*(P1-P0))
            for zb, j in JUMP.items():
                if z < zb <= z+h: T *= 1+j
            z += h
        out.append(T)
    return np.array(out)
