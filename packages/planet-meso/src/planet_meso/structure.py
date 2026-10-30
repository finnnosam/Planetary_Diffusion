"""Height encoding and conservative resampling used by the meso hand-off and its reviews.

Copied unchanged on 29 Oct 2026 from planet_training/terrain_structure.py (Planetary-Infinite-Diffusion, tag
pre-split-2026-10-29); the rest of that module belonged to the learned-layout models and was not carried over."""
from functools import lru_cache

import numpy as np
import torch
import torch.nn.functional as F


def decode(x):
    return x.sign()*x.square()


def encode(x):
    return x.sign()*x.abs().sqrt()


@lru_cache(maxsize=8)
def area_weights(size, spacing):
    """Tangent-plane cell-centre solid-angle Jacobian; rotation invariant."""
    q = (np.arange(size)+.5-size/2)*spacing/6371000.
    return torch.from_numpy((1+q[:, None]**2+q[None, :]**2)**-1.5).float()[None, None]


def overlap_matrix(old, new):
    """Exact cell interval overlaps, including noninteger 7:4 grid conversion."""
    lo = torch.arange(new, dtype=torch.float64)[:, None]*old/new
    hi = (torch.arange(new, dtype=torch.float64)[:, None]+1)*old/new
    j = torch.arange(old, dtype=torch.float64)[None]
    return (torch.minimum(hi, j+1)-torch.maximum(lo, j)).clamp_min(0).float()


def conservative_resample(z, new_size, spacing=23040.):
    """Area-average existing metre cell means; never interpolate encoded heights."""
    old = z.shape[-1]
    if z.shape[-2] != old: raise ValueError('Expected square source fields')
    a = overlap_matrix(old, new_size).to(z)
    w = area_weights(old, spacing).to(z)[0, 0]
    return (a @ (z*w) @ a.T)/(a @ w @ a.T)


def blur(x, sigma_cells):
    r = int(np.ceil(3*sigma_cells))
    q = torch.arange(-r, r+1, device=x.device, dtype=x.dtype)
    k = torch.exp(-.5*(q/sigma_cells).square()); k /= k.sum()
    x = F.conv2d(F.pad(x, (r, r, 0, 0), mode='reflect'), k[None, None, None, :])
    return F.conv2d(F.pad(x, (0, 0, r, r), mode='reflect'), k[None, None, :, None])
