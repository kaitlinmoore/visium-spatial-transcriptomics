"""Tests for the Visium loader and its coordinate-frame / scalefactor guards.

load_visium reads a real Space Ranger directory (no committed data), so it is
exercised via a monkeypatched reader; the pure seams — extract_scalefactors and
assert_spatial_frame — are tested directly. These guards are the first line
against the coordinate-frame silent-bug class.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from synthetic import make_synthetic_visium


# --- assert_spatial_frame (guard) --------------------------------------------

def test_assert_spatial_frame_passes_on_valid(adata):
    from visium_spatial.load_visium import assert_spatial_frame

    assert assert_spatial_frame(adata) is None  # synthetic fixture records the frame


def test_assert_spatial_frame_missing_coords():
    from visium_spatial.load_visium import assert_spatial_frame

    a = ad.AnnData(np.zeros((3, 1), dtype=float))
    with pytest.raises(ValueError, match="missing"):
        assert_spatial_frame(a)


def test_assert_spatial_frame_wrong_frame():
    from visium_spatial.load_visium import SPATIAL_FRAME, assert_spatial_frame

    a = ad.AnnData(np.zeros((3, 1), dtype=float))
    a.obsm["spatial"] = np.zeros((3, 2))
    a.uns["spatial_frame"] = "array_rowcol"  # wrong units recorded
    with pytest.raises(ValueError, match=SPATIAL_FRAME):
        assert_spatial_frame(a)


def test_assert_spatial_frame_wrong_shape():
    from visium_spatial.load_visium import SPATIAL_FRAME, assert_spatial_frame

    a = ad.AnnData(np.zeros((3, 1), dtype=float))
    a.obsm["spatial"] = np.zeros((3, 3))  # not (n, 2)
    a.uns["spatial_frame"] = SPATIAL_FRAME
    with pytest.raises(ValueError, match=r"\(n, 2\)"):
        assert_spatial_frame(a)


# --- extract_scalefactors -----------------------------------------------------

def test_extract_scalefactors_reads_squidpy_layout(scalefactors):
    from visium_spatial.load_visium import extract_scalefactors

    a = ad.AnnData(np.zeros((2, 1), dtype=float))
    a.uns["spatial"] = {"libX": {"scalefactors": scalefactors}}
    out = extract_scalefactors(a)
    assert out["tissue_hires_scalef"] == pytest.approx(0.15)
    assert all(isinstance(v, float) for v in out.values())


def test_extract_scalefactors_requires_single_library(scalefactors):
    from visium_spatial.load_visium import extract_scalefactors

    a = ad.AnnData(np.zeros((2, 1), dtype=float))
    a.uns["spatial"] = {"libA": {"scalefactors": scalefactors}, "libB": {"scalefactors": scalefactors}}
    with pytest.raises(ValueError, match="one Visium library"):
        extract_scalefactors(a)


def test_extract_scalefactors_missing_spatial():
    from visium_spatial.load_visium import extract_scalefactors

    with pytest.raises(ValueError, match="missing"):
        extract_scalefactors(ad.AnnData(np.zeros((2, 1), dtype=float)))


# --- load_visium (monkeypatched reader) --------------------------------------

def test_load_visium_records_frame(monkeypatch, tmp_path):
    """load_visium must stamp uns['spatial_frame'] and pass the guard, without a
    real Space Ranger directory."""
    import squidpy as sq

    from visium_spatial.load_visium import SPATIAL_FRAME, load_visium

    monkeypatch.setattr(sq.read, "visium", lambda path, **kw: make_synthetic_visium())
    adata = load_visium(tmp_path)  # tmp_path exists, so the existence check passes
    assert adata.uns["spatial_frame"] == SPATIAL_FRAME
    assert "spatial" in adata.obsm


def test_load_visium_missing_dir_raises():
    from visium_spatial.load_visium import load_visium

    with pytest.raises(FileNotFoundError):
        load_visium("does/not/exist")
