"""The plate layer's public names in one place (the code lives in the modules of planet_tectonics):

  params         PARAMS (every tunable, with source and history), boundary codes, Earth's stable-crust quantiles
  sim            PlateSim; its methods by mechanism in start, composition, crust, surface, motion, events
  plate          Plate (one rigid plate and its raster)
  isostasy, graphops, crust_kernels, boundaries   the numerical kernels

  from planet_tectonics.plates import PlateSim, PARAMS
"""
from planet_grid import fields  # re-exported: scripts use plates.fields

from planet_tectonics.params import PARAMS, RIDGE, SUBDUCTION, COLLISION, TRANSFORM, CONVERGENT, EARTH_STABLE_CRUST_Q, EARTH_STABLE_CRUST_KM
from planet_tectonics.isostasy import gdh1_depth, bulk_thickness, sediment_column, isostatic_elevation, sea_level_for_volume
from planet_tectonics.graphops import nearest_label, components, plate_holes
from planet_tectonics.crust_kernels import quantile_map, margin_shortening, edge_spread, _zone_dijkstra, rift_pure_shear, _pure_shear_apply, spread_to_band, trench_flux_consumed
from planet_tectonics.boundaries import classify_boundaries, classify_boundaries_ring
from planet_tectonics.plate import Plate
from planet_tectonics.sim import PlateSim

__all__ = ["fields", "PARAMS", "RIDGE", "SUBDUCTION", "COLLISION", "TRANSFORM", "CONVERGENT", "EARTH_STABLE_CRUST_Q", "EARTH_STABLE_CRUST_KM", "gdh1_depth", "bulk_thickness", "sediment_column", "isostatic_elevation", "sea_level_for_volume", "nearest_label", "components", "plate_holes", "quantile_map", "margin_shortening", "edge_spread", "_zone_dijkstra", "rift_pure_shear", "_pure_shear_apply", "spread_to_band", "trench_flux_consumed", "classify_boundaries", "classify_boundaries_ring", "Plate", "PlateSim"]
