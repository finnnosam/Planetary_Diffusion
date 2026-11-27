"""Independent derivation of the radial Stokes equations for one spherical-harmonic degree (pass 48 kernel check).

Axisymmetric (m = 0 suffices: radial viscosity makes the response independent of m).
u_r = U(r) P(th), u_th = V(r) P'(th), p = Pi(r) P(th), eta = eta(r); P obeys Legendre's equation
P'' = -cot(th) P' - L P with L = l(l+1). Body force f_r = -g drho(r) P(th) (a dense anomaly pulls toward the centre).
Prints the three radial equations (continuity, r-momentum / P, th-momentum / P') and the traction components.
"""
import sympy as sp

r, th, L = sp.symbols('r theta L', positive=True)
g = sp.symbols('g', positive=True)
U, V, Pi, eta, drho = [sp.Function(n)(r) for n in ('U', 'V', 'Pi', 'eta', 'drho')]
P = sp.Function('P')(th)

ur = U*P
ut = V*sp.diff(P, th)
p = Pi*P

err = sp.diff(ur, r)
ett = sp.diff(ut, th)/r+ur/r
epp = ur/r+ut*sp.cot(th)/r
ert = sp.Rational(1, 2)*(r*sp.diff(ut/r, r)+sp.diff(ur, th)/r)

srr = -p+2*eta*err
stt = -p+2*eta*ett
spp = -p+2*eta*epp
srt = 2*eta*ert

div = sp.diff(r**2*ur, r)/r**2+sp.diff(sp.sin(th)*ut, th)/(r*sp.sin(th))
fr = sp.diff(srr, r)+sp.diff(srt, th)/r+(2*srr-stt-spp+srt*sp.cot(th))/r-g*drho*P
ft = sp.diff(srt, r)+sp.diff(stt, th)/r+((stt-spp)*sp.cot(th)+3*srt)/r

# Legendre reduction: replace P'' and P''' using P'' = -cot P' - L P
P1 = sp.Symbol('P1'); P0 = sp.Symbol('P0')
P2 = -sp.cot(th)*P1-L*P0
P3 = sp.diff(-sp.cot(th)*sp.Function('Q1')(th)-L*sp.Function('Q0')(th), th)
P3 = P3.subs({sp.diff(sp.Function('Q1')(th), th): P2, sp.diff(sp.Function('Q0')(th), th): P1}).subs(
    {sp.Function('Q1')(th): P1, sp.Function('Q0')(th): P0})


def reduce(e):
    e = e.subs(sp.diff(P, th, 3), P3).subs(sp.diff(P, th, 2), P2).subs(sp.diff(P, th), P1).subs(P, P0)
    return sp.expand(sp.simplify(sp.expand(e)))


c = reduce(div)
m_r = reduce(fr)
m_t = reduce(ft)
t_rr = reduce(srr)
t_rt = reduce(srt)
# each must be a multiple of a single angular function
cc = sp.simplify(c/P0); mr = sp.simplify(m_r/P0); mt = sp.simplify(m_t/P1)
print('continuity / P:', sp.simplify(cc), '   free of theta:', not cc.has(th) and not cc.has(P1))
print('r-mom / P     :', sp.collect(sp.expand(mr), [sp.diff(U, r, 2), sp.diff(U, r), U, sp.diff(V, r), V, Pi]),
      '   free of theta:', not mr.has(th) and not mr.has(P1))
print('th-mom / P\'   :', sp.collect(sp.expand(mt), [sp.diff(V, r, 2), sp.diff(V, r), V, sp.diff(U, r), U, Pi]),
      '   free of theta:', not mt.has(th) and not mt.has(P0))
print('sigma_rr / P  :', sp.simplify(t_rr/P0))
print('sigma_rt / P\' :', sp.simplify(t_rt/P1))
