"""CRUST1.0 (Laske et al. 2013; data/crust1): its folder, the crust-type codes it counts as oceanic, and the 1 degree
cell each graph node lies in."""
import numpy as np

from planet_grid.paths import DATA

C1 = DATA/'crust1'
OCEANIC = {'A1', 'A0', 'B-', 'V1', 'Y3'}


def crust1_cell(g):
    """Index (0..64799) of the CRUST1.0 1 degree cell every node lies in (the same mapping as earth_nodes)."""
    lat = np.degrees(np.arcsin(np.clip(g.xyz[:, 2], -1, 1))); lon = np.degrees(np.arctan2(g.xyz[:, 1], g.xyz[:, 0]))
    r = np.clip((89.5-lat+.5).astype(int), 0, 179); c = np.clip((lon+180).astype(int), 0, 359)
    return r*360+c
