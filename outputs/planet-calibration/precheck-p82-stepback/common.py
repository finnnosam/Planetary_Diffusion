"""Pass 82 step-back (journal 106): shared loaders for the land-height decomposition probes. Run from the repo root.
Planet: a snapshot's elevation split into the model's own local terms (surface.elevation in the Airy limit):
continental  z_datum = airy(T) + dens(T, sed) + sedt(sed) + cterm(lith age) + dyn + rest   (rest: flexure, water load)
oceanic      z_datum = gdh(age) + crust(T) + sedt(sed) + dyn + rest
Earth on the same n = 217 nodes: ETOPO, CRUST1.0 (thickness, sediment, type), Hoggard 2020 SL2013sv L, Seton 2020
ages, GlobSed; the same law applied to Earth's inputs leaves Earth's residual."""
import json, os
from pathlib import Path
import numpy as np
from planet_grid.paths import DATA, OUTPUTS
from planet_grid.sphere_graph import CubeSphere
from planet_tectonics import lithosphere as lt
from planet_tectonics.isostasy import gdh1_depth, sediment_column, bulk_thickness, sea_level_for_volume
try:                                     # pass 82: the floor's law is the run's (ocean_depth_rule)
    from planet_tectonics.isostasy import ocean_depth
except ImportError:
    ocean_depth = lambda age, p: gdh1_depth(age)
from planet_tectonics.deep_crust import densified_mass
from planet_reviews.crust_review import run_params

ARM = 'outputs/planet-calibration/sweep-20261109-222925/lith_frame-lamb'
SEEDS = range(1001, 1007)
CACHE = Path(os.environ['P82_CACHE']) if os.environ.get('P82_CACHE') else None   # optional cache folder for earth()
_g = None
def graph():
    global _g
    if _g is None: _g = CubeSphere(217)
    return _g

def wq(x, w, q):
    """Area-weighted quantile(s)."""
    x = np.asarray(x, float); w = np.asarray(w, float)
    o = np.argsort(x); c = np.cumsum(w[o]); c = (c-.5*w[o])/c[-1]
    return np.interp(q, c, x[o])
def wmean(x, w): return float(np.sum(x*w)/np.sum(w))

def planet(seed, tag='t+0000.0', arm=ARM):
    g = graph(); d = np.load(f'{arm}/seed{seed}/snapshots/{tag}.npz'); p = run_params(f'{arm}/seed{seed}')
    rc, rm, rw = p['rho_crust'], p['rho_mantle'], p['rho_water']
    sl = float(d['sea_level_m']); zd = d['elevation_m'].astype(float)+sl
    cont = d['continental'].astype(bool); T = d['crust_thickness_m'].astype(float); sed = d['sediment_m'].astype(float)
    sh, sm_dry, sm_wet = sediment_column(sed, cont, p)
    age = d['lithosphere_age_myr'].astype(float); oa = d['ocean_age_myr'].astype(float)
    dyn = d['dynamic_topography_m'].astype(float)
    airy = p['cont_ref_elev_m']+(T-p['cont_ref_thick_m'])*(rm-rc)/rm
    dens = -densified_mass(T, np.where(cont, sm_dry, 0.), cont, p)/rm
    sedc = sh-sm_dry/rm
    cterm = np.where(cont, lt.elevation_term(age, p), 0.)
    gdh = -ocean_depth(oa, p); crust_o = (T-p['ocean_thick_m'])*(rm-rc)/(rm-rw); sedo = sh-sm_wet/(rm-rw)
    local = np.where(cont, airy+dens+sedc+cterm, gdh+crust_o+sedo)+dyn
    L = np.where(cont, lt.thickness_km(np.nan_to_num(age), p), np.nan)
    return dict(g=g, a=g.area, p=p, sea=sl, zd=zd, cont=cont, T=T, sed=sed, sh=sh, age=age, oa=oa, dyn=dyn, airy=airy,
                dens=dens, sedc=sedc, cterm=cterm, gdh=gdh, crust_o=crust_o, sedo=sedo, local=local, rest=zd-local,
                L=L, tot=(T+sh)/1e3, orog=d['orogeny_age_myr'].astype(float), tect=d['tectonic_thickening_m'].astype(float),
                btype=d['boundary_type'], dyn_air=d['dynamic_topography_air_m'].astype(float))

def earth(p=None):
    """Earth on the nodes (cached): z (ETOPO), ice, CRUST1.0 total / sediment / type, continental domain (non-oceanic
    type, total >= 20 km), Hoggard L, class age, Seton age, GlobSed; and the model's law on these inputs."""
    g = graph(); f = CACHE/'earth.npz' if CACHE else None
    if f is not None and f.exists():
        e = dict(np.load(f, allow_pickle=False))
    else:
        from planet_earth.reference import ice_mask
        lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
        r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
        c1 = DATA/'crust1'
        T = np.loadtxt(c1/'crsthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
        S = np.loadtxt(c1/'sedthk.xyz')[:, 2].reshape(180, 360)[r, c]*1e3
        tn = np.array([l.split() for l in (c1/'CNtype1-1.txt').read_text().splitlines() if l.strip()])[r, c]
        dd = np.loadtxt(DATA.joinpath(*lt.HOGGARD_SL2013SV), delimiter=',')
        grid = np.full((361, 721), np.nan)
        grid[np.round((90-dd[:, 1])*2).astype(int), np.round((dd[:, 0] % 360)*2).astype(int) % 720] = dd[:, 2]
        grid[:, 720] = grid[:, 0]
        L = grid[np.round((90-lat)*2).astype(int), np.round((lon % 360)*2).astype(int)]
        R = OUTPUTS/'planet-reference'
        e = dict(z=np.load(R/'earth-n217.npz')['elevation_m'].astype(float), ice=ice_mask(g), T=T, S=S, tn=tn.astype('U2'),
                 L=L, cage=np.array([lt.CLASS_AGE_MA.get(t, np.nan) for t in tn]), oceanic=np.isin(tn, lt.OCEANIC_TYPES),
                 oa=np.load(R/'earth-seafloor-age-n217.npz')['age_myr'].astype(float),
                 gsed=np.load(R/'globsed-n217.npz')['sediment_m'].astype(float), lat=lat, lon=lon)
        if f is not None:
            CACHE.mkdir(parents=True, exist_ok=True); np.savez(f, **e)
    if p is None:
        p = run_params(f'{ARM}/seed1001')
    rc, rm, rw, rs = p['rho_crust'], p['rho_mantle'], p['rho_water'], p['rho_sediment']
    e['g'] = g; e['a'] = g.area; e['p'] = p
    e['cont'] = ~e['oceanic'] & (e['T'] >= 20e3)
    Tc = e['T']-e['S']; e['Tc'] = Tc; e['tot'] = e['T']/1e3
    e['airy'] = p['cont_ref_elev_m']+(Tc-p['cont_ref_thick_m'])*(rm-rc)/rm
    # sediment as the model books it: CRUST1.0's thickness is the bulk height of the model's compacting column
    # (continental law), its dry mass the load; the uniform rho_sediment reading (coast_budget's) kept beside it
    h = e['S']; phi0, c = p['phi0_continent'], p['compaction_c_continent_per_km']/1e3
    solid = h-phi0/c*(1-np.exp(-c*h)); sm_dry = p['rho_grain']*solid+rw*(h-solid)
    e['dens'] = -densified_mass(Tc, sm_dry, np.ones(Tc.size, bool), p)/rm
    e['sedc'] = h-sm_dry/rm
    e['sedc_2400'] = h*(1-rs/rm)
    e['c_age'] = lt.elevation_term(e['cage'], p)
    k = lt.LAMB_K_M_PER_KM
    e['c_seis'] = -k*(e['L']-lt.LAMB_REF_KM)
    e['c_seis_clip'] = -k*(np.clip(e['L'], lt.LAMB_L_KM[0], lt.LAMB_L_KM[-1])-lt.LAMB_REF_KM)
    e['zd'] = e['z']; e['sea'] = 0.
    return e
