"""Isostasy and sea level: ocean depth from age (GDH1), sediment compaction, Airy elevation, the sea level that holds
the water inventory."""
import numpy as np

from planet_tectonics.params import PARAMS



def gdh1_depth(age_myr):
    """Ocean depth (m, positive) from lithosphere age (Stein & Stein 1992)."""
    a = np.clip(age_myr, 0, None)
    return np.where(a < 20, 2600.+365.*np.sqrt(a), 5651.-2473.*np.exp(-.0278*a))



def bulk_thickness(S, phi0, c_per_km):
    """Bulk thickness h (m) of a sediment column holding solid thickness S (m) when its porosity falls with depth z
    below the top as phi0 exp(-c z) (Athy 1930): S = h - phi0 / c (1 - exp(-c h)). Newton's method from the thin-cover
    bound S / (1 - phi0), from which it converges monotonically (the residual is convex and increasing in h)."""
    S = np.maximum(np.asarray(S, np.float64), 0.)
    phi0 = np.asarray(phi0, np.float64); c = np.asarray(c_per_km, np.float64)/1e3
    h = S/(1-phi0)
    for _ in range(50):
        e = np.exp(-c*h)
        f = h-phi0/c*(1-e)-S
        h = h-f/(1-phi0*e)
        if np.all(np.abs(f) < 1e-6):
            break
    return h



def sediment_column(sed, cont, p=PARAMS):
    """Sediment column of every node: (surface height, mass per area dry, mass per area submerged minus the water it
    displaces). With sediment_compaction, sed is solid thickness: bulk height from bulk_thickness (shale law on oceanic
    crust, shaly sandstone on continental crust), mass rho_grain S + rho_water (h - S) with water-filled pores, and
    (rho_grain - rho_water) S under water. Without it, a layer of rho_sediment: (sed, rho_sediment sed,
    (rho_sediment - rho_water) sed)."""
    rw = p['rho_water']
    if not p.get('sediment_compaction'):
        rs = p['rho_sediment']
        return sed, rs*sed, (rs-rw)*sed
    phi0 = np.where(cont, p['phi0_continent'], p['phi0_ocean'])
    c = np.where(cont, p['compaction_c_continent_per_km'], p['compaction_c_ocean_per_km'])
    h = bulk_thickness(sed, phi0, c)
    rg = p['rho_grain']
    return h, rg*sed+rw*(h-sed), (rg-rw)*sed



def isostatic_elevation(cont, thick, ocean_age, p=PARAMS):
    """Elevation (m): oceanic crust from age; continental crust from Airy thickness, water-loaded below 0 m."""
    dry = p['cont_ref_elev_m']+(thick-p['cont_ref_thick_m'])*(p['rho_mantle']-p['rho_crust'])/p['rho_mantle']
    wet = np.where(dry < 0, dry*p['rho_mantle']/(p['rho_mantle']-p['rho_water']), dry)
    return np.where(cont, wet, -gdh1_depth(ocean_age))



def sea_level_for_volume(z, area, volume_m3):
    """Level s (m) at which the water below it, sum of (s - z) * area over z < s, equals volume_m3."""
    order = np.argsort(z)
    zs, a = z[order], area[order]
    ca = np.cumsum(a); cza = np.cumsum(zs*a)
    held = zs*ca-cza                                  # volume below level zs[k] (the k+1 lowest nodes flooded)
    k = int(np.searchsorted(held, volume_m3))
    if k == 0:
        return float(zs[0])
    k = min(k, len(zs))
    return float((volume_m3+cza[k-1])/ca[k-1])
