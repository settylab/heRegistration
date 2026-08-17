"""Test bootstrap: put src/ on sys.path so `import hexenium` resolves
without a `pip install -e .`. Lets smoke tests run in a stock Python
env with just pytest + pyyaml installed."""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
