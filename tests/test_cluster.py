"""Tests for Leiden clustering in expression space (fresh AnnData per test).

These assert the mechanics (every spot labeled, key_added honored), not a
specific cluster count — cluster granularity on a 3-gene toy is not meaningful.
The important conceptual point (expression graph, not spatial) is enforced by
the module using sc.pp.neighbors, verified here only in that it does not touch
obsp["spatial_connectivities"].
"""

from __future__ import annotations

import pytest

from synthetic import make_synthetic_visium


def test_leiden_labels_every_spot():
    from visium_spatial.cluster import leiden_clusters
    from visium_spatial.preprocess import normalize

    a = make_synthetic_visium(seed=0)
    normalize(a)
    leiden_clusters(a, resolution=1.0, n_neighbors=15)

    assert "leiden" in a.obs
    assert a.obs["leiden"].notna().all()
    assert a.obs["leiden"].nunique() >= 1


def test_leiden_key_added_and_not_spatial_graph():
    from visium_spatial.cluster import leiden_clusters
    from visium_spatial.preprocess import normalize

    a = make_synthetic_visium(seed=0)
    normalize(a)
    leiden_clusters(a, key_added="compartment")

    assert "compartment" in a.obs
    # expression clustering must not have built the spatial weights graph
    assert "spatial_connectivities" not in a.obsp
