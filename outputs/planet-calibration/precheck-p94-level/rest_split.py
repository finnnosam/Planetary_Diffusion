"""Pass 94, probe 2: what the planet's 'rest' (z - the local column terms) holds, by crust bin, on today's control.
surface.elevation adds the water load of a thinned crust's own Airy depth, own_airy * rw / (rm - rw), wherever
own_airy < 0 (line 180), also where the column's sediment reaches above sea level and is booked as a dry load
(c_dry): there no water fills the column. Split of rest per bin (area means, six seeds, 0 Myr):
  w_crust   the booked crust water term (all own_airy < 0)
  w_dry     its part on c_dry columns (no water above them in the model's own load)
  w_C       C's water amplification (C * rw / (rm - rw) where the surface is under the sea)
  flex      rest - w_crust - w_C (flexure: the deflection's departure from local Airy, plus anything else)
and the local column with the 1-D Airy water law applied once, alike for planet and Earth: z_air = airy + dens + sedc
+ C (+ dyn on the planet); z_loc = z_air * rm / (rm - rw) where z_air < sea, else z_air. Earth's residual z - z_loc
beside the planet's z - z_loc. Earth's dyn left in its residual unless h5py reads Holdt's field. The booked term follows the run's water_load_rule."""
import sys, numpy as np
sys.path.insert(0, 'outputs/planet-calibration/precheck-p82-stepback')
import common as C
from planet_tectonics.isostasy import sediment_column
ARM = sys.argv[1] if len(sys.argv) > 1 else 'outputs/planet-calibration/sweep-20261119-182758/backarc_extension-True_rollback_rule-fixed_rollback_mm_yr-9'
B = [(20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 55), (30, 45), (20, 80)]
e = C.earth(C.run_params(f'{ARM}/seed1001')); p = e['p']; rm, rw = p['rho_mantle'], p['rho_water']
try:
    import h5py
    Dh = 'data/residual_topography/Residual depthtopographic measurements/Residual depthtopographic measurements/Holdt_et_al_2022_SupportingInformation/'
    with h5py.File(Dh+'degree_1_to_30.grd', 'r') as f: lon, lat, v = f['lon'][:], f['lat'][:], f['z'][:]
    i = np.clip(np.round((e['lat']-lat[0])/(lat[1]-lat[0])).astype(int), 0, lat.size-1)
    j = np.round(((e['lon'] % 360)-lon[0])/(lon[1]-lon[0])).astype(int) % lon.size
    edyn = np.nan_to_num(v[i, j]*1e3); dyn_note = 'Holdt 1-30 taken out'
except Exception as x:
    edyn = np.zeros(e['z'].size); dyn_note = f'dyn left in ({x})'
def loc(zair, sea): return np.where(zair < sea, sea+(zair-sea)*rm/(rm-rw), zair)
ez_air = e['airy']+e['dens']+e['sedc']+e['c_seis_clip']+edyn; ez_loc = loc(ez_air, 0.)
dom_e = e['cont'] & ~e['ice'] & np.isfinite(e['L'])
rows = {b: [] for b in B}; dryshare = {b: [] for b in B}
for s in C.SEEDS:
    P = C.planet(s, 't+0000.0', ARM); a = P['a']; sl = P['sea']
    own = P['airy']+P['dens']-sl
    if P['p'].get('water_load_rule', 'parts') == 'column':     # the run's own rule (pass 94)
        smd = sediment_column(P['sed'], P['cont'], P['p'])[1]
        wet = P['cont'] & (own+P['sh']-np.where(P['cont'], smd, 0.)/rm < 0); c_dry = P['cont'] & ~wet
    else:
        c_dry = P['cont'] & (np.where(own < 0, own*rm/(rm-rw), own)+P['sh'] > 0); wet = P['cont']
    w_crust = np.where(wet & (own < 0), own*rw/(rm-rw), 0.)
    w_dry = np.where(c_dry, w_crust, 0.)
    w_C = np.where(P['cont'] & (P['zd'] < sl), P['cterm']*rw/(rm-rw), 0.)
    flex = P['rest']-w_crust-w_C
    zair = P['airy']+P['dens']+P['sedc']+P['cterm']+P['dyn']; zloc = loc(zair, sl)
    for b in B:
        m = P['cont'] & (P['tot'] >= b[0]) & (P['tot'] < b[1]); am = a[m]
        rows[b].append([C.wmean(x[m], am) for x in (P['rest'], w_crust, w_dry, w_C, flex, P['zd']-zloc, zloc, P['zd'])])
        dryshare[b].append(am[c_dry[m] & (own[m] < 0)].sum()/am.sum())
print(ARM, '|', dyn_note)
print('bin km | planet: rest = w_crust (of which w_dry) + w_C + flex | share of the bin dry with own_airy < 0 | planet z - z_loc | Earth z - z_loc | z_loc planet / Earth | z planet / Earth (m)')
for b in B:
    r = np.array(rows[b]); mm = r.mean(0); sd = r.std(0, ddof=1)
    me = dom_e & (e['tot'] >= b[0]) & (e['tot'] < b[1]); ae = e['a'][me]
    print(f'  {b[0]:2d}-{b[1]:2d} | {mm[0]:+5.0f} = {mm[1]:+5.0f} ({mm[2]:+5.0f}) + {mm[3]:+4.0f} + {mm[4]:+5.0f} | {np.mean(dryshare[b]):.2f} | {mm[5]:+5.0f} ({sd[5]:3.0f}) | {C.wmean((e["z"]-ez_loc)[me], ae):+5.0f} | {mm[6]:+5.0f} / {C.wmean(ez_loc[me], ae):+5.0f} | {mm[7]:+5.0f} / {C.wmean(e["z"][me], ae):+5.0f}')
