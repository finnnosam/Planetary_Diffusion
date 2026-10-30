"""Lets pytest import the stage packages from packages/*/src without installing them (an editable install via
requirements-dev.txt does the same for scripts)."""
import sys
from pathlib import Path

for src in sorted(Path(__file__).parent.glob('packages/*/src')):
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))
