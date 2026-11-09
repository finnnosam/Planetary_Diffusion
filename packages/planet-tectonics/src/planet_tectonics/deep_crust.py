"""Deep-crust densification with pressure (seventy-sixth pass, journal 100): `crust_density_rule` 'pressure'.

The model's crust has one density (rho_crust 2,830 kg/m^3, CRUST1.0's crystalline mean 2,827), so a thickened column
stays as buoyant per km at 60 km as at 35 km. Earth's lower crust is mafic granulite (CRUST1.0's crust types: the
lowest layer is 0.324 of the crystalline crust at 2,947 kg/m^3, area means over the continental types; the share is
the same in every thickness bin, 30-55 km). Pushed deeper by thickening, mafic granulite turns to garnet granulite and
then eclogite: garnet comes in at ~1.0 GPa (gabbro -> garnet granulite, Green & Ringwood 1967, experimental) and the
rock reaches the mantle's density at ~2 GPa (Chapman et al. 2019; the same source as `eclogite_rule`'s 65 km, where
the root founders). Here the lower crust (the deepest lower_crust_share of the crystalline column) gains density
linearly in pressure from rho_lower_crust at densify_start_gpa to rho_mantle at densify_neutral_gpa, and holds it
beyond (crust past eclogite_thick_m founders under `eclogite_rule`). Equilibrium is assumed: dry granulite can stay
metastable (Jackson et al. 2004), so this is the upper bound of the effect.

Pressure at depth z below the crust's top: g (sediment dry mass per area + rho_crust z). The excess mass per area,
M = integral over the lower crust of drho(P(z)) dz (kg/m^2), is exact for the linear ramp: with u(z) the ramp's
fraction, M = drho_max (G(u(base)) - G(u(top of lower crust))) / (du/dz), G(u) = 0 (u <= 0), u^2/2 (0 < u <= 1),
u - 1/2 (u > 1). surface.elevation adds M to the continental crust's dry flexural load (long wavelengths: the surface
sinks M / rho_mantle); elevation_check's own Airy level subtracts M / rho_mantle. Volume and thickness are unchanged:
the crust's ledgers do not see it.
"""
import numpy as np

G_EARTH = 9.81


def _ramp_integral(u):
    return np.where(u <= 0., 0., np.where(u <= 1., .5*u*u, u-.5))


def densified_mass(thick, sed_dry_mass, cont, p):
    """Excess mass per area (kg/m^2) of the densified lower crust; 0 off the continent and under 'uniform'."""
    thick = np.asarray(thick, np.float64)
    if p.get('crust_density_rule', 'uniform') == 'uniform':
        return np.zeros_like(thick)
    if p['crust_density_rule'] != 'pressure':
        raise ValueError(f"crust_density_rule {p['crust_density_rule']!r}: 'uniform' or 'pressure'")
    rc = p['rho_crust']
    p1, p2 = p['densify_start_gpa']*1e9, p['densify_neutral_gpa']*1e9
    drho = p['rho_mantle']-p['rho_lower_crust']
    over = np.maximum(np.asarray(sed_dry_mass, np.float64), 0.)
    T = np.maximum(thick, 0.)
    k = G_EARTH*rc/(p2-p1)                                  # du/dz (1/m)
    u = lambda z: (G_EARTH*(over+rc*z)-p1)/(p2-p1)
    z0 = (1.-p['lower_crust_share'])*T
    m = drho*(_ramp_integral(u(T))-_ramp_integral(u(z0)))/k
    return np.where(np.asarray(cont, bool), m, 0.)
