"""Where the Earth datasets and the results live. Neither is part of the repository.

PLANET_DATA (default: ./data) holds the downloaded Earth datasets; PLANET_OUTPUTS (default: ./outputs) holds runs,
sweeps, reviews and model checkpoints. Relative values resolve against the current folder when this module is first
imported, so run commands from the repository root or set the variables."""
import os
from pathlib import Path

DATA = Path(os.environ.get('PLANET_DATA', 'data')).resolve()
OUTPUTS = Path(os.environ.get('PLANET_OUTPUTS', 'outputs')).resolve()
