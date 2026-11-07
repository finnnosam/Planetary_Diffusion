"""K (sixty-fifth pass, journal 87): Earth's rifted-margin production since its supercontinent, for split_rule 'earth'.

EARTH_RIFT_PRODUCTION rows (t Ma, km of rifted margin made per Myr over (t, t+10] Ma, Earth's non-ocean area Mkm^2),
measured on the n = 217 graph from Mueller et al. (2019): at each 10 Myr reconstruction, the continental nodes (no
paleo-age = not ocean) beside ocean floor <= 10 Myr old, leaving out those within 300 km of a subduction zone resolved
from the same plate model with pygplates (their young floor is a back-arc or a young plate at a trench); length =
nodes x graph spacing; floor that young beside a margin was made in the last 10 Myr, so each margin is counted at its
birth, both conjugate sides (outputs/planet-calibration/precheck-p65-towers/earth_rift_production.py). Mean 1,553
km/Myr (7.4 per Mkm^2), 388k km in all; low while Pangea was whole (250-170 Ma), the breakup peak at 160-110 Ma.
The planet's start (a supercontinent) is mapped onto 250 Ma; production is per Mkm^2 of continent, so a planet
with less continent makes proportionally less margin. Runs longer than 250 Myr keep the oldest row before the start
and the youngest after 0 Ma.
"""
import numpy as np

EARTH_RIFT_PRODUCTION = (
    (0, 1423, 218.7), (10, 1309, 218.2), (20, 1763, 217.6), (30, 1542, 217.7), (40, 1045, 219.1),
    (50, 2226, 217.9), (60, 1997, 216.3), (70, 1517, 214.9), (80, 2915, 214.0), (90, 680, 211.8),
    (100, 948, 210.5), (110, 2817, 210.2), (120, 3463, 208.7), (130, 2694, 206.7), (140, 3034, 203.8),
    (150, 4308, 203.4), (160, 2311, 202.0), (170, 208, 201.4), (180, 446, 200.2), (190, 569, 201.7),
    (200, 0, 201.3), (210, 0, 199.5), (220, 242, 199.8), (230, 399, 201.6), (240, 977, 197.4))
START_MA = 250.
_PER_MKM2 = np.array([r[1]/r[2] for r in EARTH_RIFT_PRODUCTION])     # km / Myr / Mkm^2, index t Ma / 10


def production_per_mkm2(t_since_start_myr):
    """Earth's rifted-margin production (km of conjugate margin per Myr per Mkm^2 of continent) at the time since the
    start: Earth's age START_MA - t_since_start, in its 10 Myr bin ((t, t+10] Ma -> row t)."""
    ma = START_MA-float(t_since_start_myr)
    i = int(np.clip(np.ceil(ma/10.-1e-9)-1, 0, len(_PER_MKM2)-1))
    return float(_PER_MKM2[i])
