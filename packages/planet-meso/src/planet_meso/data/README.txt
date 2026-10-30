synthetic_map_stats.json and procedural.py (one folder up) were copied unchanged on 29 Oct 2026 from the original
generator (github.com/finnnosam/Planetary-Infinite-Diffusion, tag pre-split-2026-10-29, planet_diffusion/).

synthetic_map_stats.json was computed by Terrain Diffusion's inference/synthetic_map.py (github.com/xandergos/
terrain-diffusion, commit e8dcb4b1a834ab2f6b1a6f5256ed7c9f2f3e8230) with frequency_mult=[1.5,3,3,3,3] and
drop_water_pct=0.5, from its upstream/data/global/etopo_10m.tif and WorldClim 2.1 10-arc-minute BIO1, BIO4, BIO12,
BIO15. Terrain Diffusion's licence (MIT) is in terrain_diffusion_LICENSE.txt. The meso stage uses its 2-D noise
quantiles for the random-terrain input source.
