"""The Earth datasets the planet is compared with: what each is, where to get it, which files the code reads, and
their SHA-256 checksums. The data are not part of the repository; put them under data/ (or set PLANET_DATA).

  python -m planet_earth.manifest            list the datasets and check that the files are there (and their checksums)
  python -m planet_earth.manifest --record   compute the checksums of the files present and save them with the package
"""
import argparse
import hashlib
import json
from pathlib import Path

from planet_grid.paths import DATA

CHECKSUMS = Path(__file__).with_name('data_checksums.json')

DATASETS = [
    {'key': 'etopo', 'name': 'ETOPO 2022 global relief model, 30 arc-second surface elevation (NOAA NCEI)',
     'files': ['etopo_2022_30s_surface.tif'],
     'source': 'https://www.ncei.noaa.gov/products/etopo-global-relief-model', 'doi': '10.25921/fd45-gt74',
     'used_by': 'Earth on the cube grids (planet_earth.reference): hypsometry, texture, elevation tests'},
    {'key': 'etopo_10m', 'name': 'ETOPO 2022 at 10 arc-minutes (2160 x 1080, WGS 84 + EGM2008 heights)',
     'files': ['etopo_10m.tif'],
     'source': 'a 10 arc-minute version of the ETOPO 2022 surface grid above, as used by Terrain Diffusion '
               '(https://github.com/xandergos/terrain-diffusion); the checksum pins the copy in use',
     'doi': '10.25921/fd45-gt74',
     'used_by': 'the meso training set (planet_meso.build_dataset --dem data/etopo_10m.tif)'},
    {'key': 'crust1', 'name': 'CRUST1.0 global crustal model, 1 degree (Laske et al. 2013)',
     'files': ['crust1/crsthk.xyz', 'crust1/sedthk.xyz', 'crust1/CNtype1-1.txt', 'crust1/CNtype1_key.txt'],
     'source': 'https://igppweb.ucsd.edu/~gabi/crust1.html',
     'used_by': 'continental crust thickness, margins, crust types (planet_reviews.crust_review, margin_width, interior_probe)'},
    {'key': 'crust1_summary', 'name': 'CRUST1.0 continental summary (derived on 23 Oct 2026, not a download)',
     'files': ['crust1/crust1-continental-summary.json'], 'derived': True,
     'source': 'computed from CRUST1.0 crsthk.xyz and CNtype1-1.txt (continental types, cos-latitude weights)',
     'used_by': 'planet_reviews.crust_review (Earth line), tests/test_crust_reviews.py'},
    {'key': 'bird2003', 'name': 'PB2002 plate boundaries, orogens and boundary steps (Bird 2003, G-cubed)',
     'files': ['tectonicplates-master/GeoJSON/PB2002_boundaries.json', 'tectonicplates-master/GeoJSON/PB2002_orogens.json',
               'tectonicplates-master/PB2002_steps.dbf'],
     'source': 'https://github.com/fraxen/tectonicplates', 'doi': '10.1029/2001GC000252',
     'used_by': 'orogens and boundaries (planet_earth.bird), boundary lengths and rates (planet_reviews.trench_origin)'},
    {'key': 'globsed', 'name': 'GlobSed v3 total sediment thickness of the oceans (Straume et al. 2019)',
     'files': ['GlobSed/V3/GlobSed_package3/ArcGIS_ascii/globsed-v3.asc'],
     'source': 'https://b2find.eudat.eu/dataset/3c5b69bd-5632-530b-92ab-2dbded62588f', 'doi': '10.1029/2018GC008115',
     'used_by': 'ocean-floor sediment (planet_earth.sediment, planet_reviews.ocean_sediment_probe)'},
    {'key': 'gsrm', 'name': 'Global Strain Rate Model v2.1, average strain (Kreemer et al. 2014)',
     'files': ['gsrm/GSRM_average_strain_v2.1.txt'],
     'source': 'https://geodesy.unr.edu/cornekreemer/gsrm.htm', 'doi': '10.1002/2014GC005407',
     'used_by': "Earth's deforming crust (planet_earth.deformation)"},
    {'key': 'octopus', 'name': 'OCTOPUS v2 basin-averaged 10Be erosion rates, CRN global collection (Codilean et al. 2022)',
     'files': ['octopus/crn_int_basins.json'],
     'source': 'https://essd.copernicus.org/articles/14/3695/2022/', 'doi': '10.5194/essd-14-3695-2022',
     'note': 'WFS layer be10-denude:crn_int_basins as GeoJSON, EPSG:4326',
     'used_by': 'the erosion-law check (scripts/erosion_law_check.py)'},
    {'key': 'seton2020', 'name': 'Present-day ocean-floor age grid (Seton et al. 2020)',
     'files': ['EarthByte/Grids/age.2020.1.GTS2012.2m.grd'],
     'source': 'https://www.earthbyte.org/webdav/ftp/earthbyte/agegrid/2020/Grids/', 'doi': '10.1029/2020GC009214',
     'used_by': 'ocean-floor age and production (planet_earth.seafloor_age, planet_reviews.plate_pieces_review)'},
    {'key': 'riveratlas', 'name': 'HydroATLAS RiverATLAS v1.0 (Linke et al. 2019); optional, not used yet',
     'files': [], 'optional': True, 'folder': 'RiverATLAS_Data_v10_shp',
     'source': 'https://www.hydrosheds.org/hydroatlas', 'doi': '10.6084/m9.figshare.9890531',
     'used_by': 'drainage validation of the landscape stage (planned)'},
]


def sha256(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(block), b''):
            h.update(chunk)
    return h.hexdigest()


def check(data=DATA, verify=True):
    """{dataset key: [(file, status)]}; status is 'ok', 'missing', 'changed' (checksum differs) or 'present' (no
    checksum recorded)."""
    sums = json.loads(CHECKSUMS.read_text()) if CHECKSUMS.exists() else {}
    out = {}
    for d in DATASETS:
        rows = []
        for f in d['files']:
            p = Path(data)/f
            if not p.exists():
                rows.append((f, 'missing'))
            elif verify and f in sums:
                rows.append((f, 'ok' if sha256(p) == sums[f] else 'changed'))
            else:
                rows.append((f, 'present'))
        if d.get('folder'):
            rows.append((d['folder']+'/', 'present' if (Path(data)/d['folder']).is_dir() else 'missing (optional)'))
        out[d['key']] = rows
    return out


def record(data=DATA):
    sums = {f: sha256(Path(data)/f) for d in DATASETS for f in d['files'] if (Path(data)/f).exists()}
    CHECKSUMS.write_text(json.dumps(sums, indent=1, sort_keys=True)+'\n')
    return sums


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--record', action='store_true', help='save the checksums of the files present')
    p.add_argument('--no-verify', action='store_true', help='check presence only (skip the checksums)')
    a = p.parse_args()
    if a.record:
        sums = record()
        print(f'{len(sums)} checksums written to {CHECKSUMS}')
    print(f'data folder: {DATA}')
    res = check(verify=not a.no_verify)
    bad = 0
    for d in DATASETS:
        print(f"\n{d['name']}\n  source: {d['source']}" + (f"  (doi:{d['doi']})" if d.get('doi') else '')
              + (f"\n  note: {d['note']}" if d.get('note') else '') + f"\n  used by: {d['used_by']}")
        for f, st in res[d['key']]:
            print(f'  {st:8s} {f}')
            bad += st in ('missing', 'changed')
    print(f"\n{'all required files present' if not bad else f'{bad} required file(s) missing or changed'}")
    raise SystemExit(1 if bad else 0)


if __name__ == '__main__':
    main()
