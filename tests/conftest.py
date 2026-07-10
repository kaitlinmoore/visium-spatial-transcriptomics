"""Shared pytest fixtures built on the committed synthetic Visium data.

The synthetic AnnData (tests/fixtures/synthetic.py) is deterministic, so these
fixtures are session-scoped and safe to share across tests. ``pyproject.toml``
puts ``src`` on the path, so tests import the modules under test directly
(``import build_graph``), matching how the notebooks use them.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"
sys.path.insert(0, str(FIXTURES))  # allow `import synthetic`

from synthetic import hh_patch_indices, make_synthetic_visium  # noqa: E402


@pytest.fixture(scope="session")
def adata():
    """A fresh synthetic Visium-like AnnData (seed=0)."""
    return make_synthetic_visium(seed=0)


@pytest.fixture(scope="session")
def hh_indices():
    """Integer positions of the planted High-High patch spots."""
    return hh_patch_indices(seed=0)


@pytest.fixture(scope="session")
def scalefactors():
    """The synthetic scalefactors dict from fixtures/scalefactors.json."""
    with open(FIXTURES / "scalefactors.json") as fh:
        return json.load(fh)
