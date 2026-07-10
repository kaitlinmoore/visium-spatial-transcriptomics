"""Tests for the scalefactor-aware H&E overlay.

The owned seam is the scalefactor multiply (to_hires_coords) — the single place
full-res coords become hires-image coords. The load-bearing plot test asserts the
scatter is drawn at the SCALED coordinates, which is the real guard against the
silent full-res-vs-hires misalignment CLAUDE.md warns about.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless; no display needed for tests

import numpy as np
import pytest


def test_to_hires_coords_scales():
    from visium_spatial.overlay import to_hires_coords

    spatial = np.array([[100.0, 200.0], [0.0, 50.0]])
    out = to_hires_coords(spatial, 0.15)
    assert np.allclose(out, [[15.0, 30.0], [0.0, 7.5]])


def test_to_hires_coords_uses_fixture_factor(scalefactors):
    from visium_spatial.overlay import to_hires_coords

    f = scalefactors["tissue_hires_scalef"]
    spatial = np.array([[1000.0, 2000.0]])
    assert np.allclose(to_hires_coords(spatial, f), [[1000.0 * f, 2000.0 * f]])


def test_to_hires_coords_rejects_bad_input():
    from visium_spatial.overlay import to_hires_coords

    with pytest.raises(ValueError, match="> 0"):
        to_hires_coords(np.array([[1.0, 2.0]]), 0.0)
    with pytest.raises(ValueError, match=r"\(n, 2\)"):
        to_hires_coords(np.array([1.0, 2.0, 3.0]), 0.15)


def test_plot_hotspots_draws_at_scaled_coords(adata, scalefactors):
    """The scatter offsets must equal the scalefactor-scaled coordinates, not the
    raw full-res coords — the anti-misalignment contract."""
    from visium_spatial.overlay import plot_hotspots, to_hires_coords

    f = scalefactors["tissue_hires_scalef"]
    values = np.arange(adata.n_obs, dtype=float)  # numeric -> colorbar branch
    ax = plot_hotspots(adata, values, scalefactor=f)

    offsets = ax.collections[0].get_offsets()
    expected = to_hires_coords(adata.obsm["spatial"], f)
    assert np.allclose(np.asarray(offsets), expected)


def test_plot_hotspots_sets_image_background(adata, scalefactors):
    from visium_spatial.overlay import plot_hotspots

    image = np.zeros((80, 100, 3), dtype=float)  # dummy hires backdrop
    ax = plot_hotspots(adata, np.arange(adata.n_obs), scalefactor=scalefactors["tissue_hires_scalef"], image=image)
    assert len(ax.images) == 1


def test_plot_hotspots_accepts_categorical_labels(adata, scalefactors):
    """LISA quadrant labels (strings) route through the legend branch."""
    from visium_spatial.overlay import plot_hotspots

    labels = np.array(["HH", "ns"] * (adata.n_obs // 2) + ["ns"] * (adata.n_obs % 2))
    ax = plot_hotspots(adata, labels, scalefactor=scalefactors["tissue_hires_scalef"])
    assert ax.get_legend() is not None
    assert len(ax.collections[0].get_offsets()) == adata.n_obs


def test_plot_hotspots_rejects_length_mismatch(adata, scalefactors):
    from visium_spatial.overlay import plot_hotspots

    with pytest.raises(ValueError, match="length"):
        plot_hotspots(adata, np.arange(adata.n_obs - 1), scalefactor=scalefactors["tissue_hires_scalef"])
