"""Load a 10x Visium section and surface its coordinate frame + scalefactors.

Wraps ``squidpy.read.visium`` and pulls out the two things the rest of the
pipeline must never guess about (see the coordinate/scalefactor hygiene section
of CLAUDE.md):

- ``obsm["spatial"]`` is in **full-resolution image pixel coordinates** (x, y).
  Array (row/col) indexing and y-axis orientation differ between image and plot
  conventions, so we record which frame is in use rather than assume it.
- ``scalefactors_json.json`` carries ``tissue_hires_scalef`` /
  ``tissue_lowres_scalef``; overlay.py multiplies the full-res coordinates by the
  hires factor to land on ``tissue_hires_image.png``. Mixing full-res coordinates
  with a scaled image is the tissue-image analog of a CRS/projection mismatch.

Design choice: loading is deliberately dumb — no QC, no normalization, no graph.
Those are separate, separately testable stages so a bad default here cannot
silently corrupt the coordinate frame that everything downstream depends on.

Alternative considered: fold QC/normalization into the loader for convenience.
Rejected — it couples the one irreversible step (fixing the coordinate frame) to
mutable analysis choices, exactly the coupling the hygiene rules warn against.
"""

from __future__ import annotations

from pathlib import Path

SPATIAL_FRAME = "fullres_pixels"  # what obsm["spatial"] is in, recorded in uns


def load_visium(outs_dir: str | Path, *, count_file: str = "filtered_feature_bc_matrix.h5"):
    """Read a Space Ranger ``outs/`` directory into an AnnData.

    Asserts the spatial coordinate frame is present and records it in
    ``adata.uns["spatial_frame"]`` so downstream code can check, not assume.

    Parameters
    ----------
    outs_dir:
        Space Ranger ``outs`` directory (contains ``spatial/`` + the matrix .h5).
    count_file:
        Name of the count-matrix HDF5 to read.

    Returns
    -------
    anndata.AnnData
    """
    import squidpy as sq

    outs_dir = Path(outs_dir)
    if not outs_dir.exists():
        raise FileNotFoundError(f"Space Ranger outs dir not found: {outs_dir}")

    adata = sq.read.visium(outs_dir, counts_file=count_file)
    # Record the coordinate frame explicitly (the one irreversible fact), then
    # verify it before returning so a malformed read fails here, not downstream.
    adata.uns["spatial_frame"] = SPATIAL_FRAME
    assert_spatial_frame(adata)
    return adata


def extract_scalefactors(adata) -> dict:
    """Return the Visium scalefactors as a plain ``dict`` of floats.

    Reads squidpy's ``uns["spatial"][library_id]["scalefactors"]`` layout. Keys:
    ``spot_diameter_fullres``, ``tissue_hires_scalef``, ``tissue_lowres_scalef``,
    ``fiducial_diameter_fullres``.
    """
    spatial = adata.uns.get("spatial")
    if not spatial:
        raise ValueError('adata.uns["spatial"] missing; not a squidpy-read Visium AnnData')

    libraries = list(spatial.keys())
    if len(libraries) != 1:
        raise ValueError(f"expected exactly one Visium library, found {libraries}")

    scalefactors = spatial[libraries[0]].get("scalefactors")
    if scalefactors is None:
        raise ValueError(f'no scalefactors under uns["spatial"]["{libraries[0]}"]')
    return {key: float(value) for key, value in scalefactors.items()}


def assert_spatial_frame(adata) -> None:
    """Fail loudly if ``obsm["spatial"]`` is missing or the frame is unrecorded.

    A guard, not a transform: catches the silent-bug class where coordinates are
    absent, in the wrong units, or of the wrong shape before any statistic runs.
    """
    if "spatial" not in adata.obsm:
        raise ValueError('adata.obsm["spatial"] is missing; no spatial coordinates')

    xy = adata.obsm["spatial"]
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f'obsm["spatial"] must be (n, 2); got shape {xy.shape}')
    if xy.shape[0] != adata.n_obs:
        raise ValueError(
            f'obsm["spatial"] has {xy.shape[0]} rows but n_obs is {adata.n_obs}'
        )

    frame = adata.uns.get("spatial_frame")
    if frame != SPATIAL_FRAME:
        raise ValueError(f'uns["spatial_frame"] must be {SPATIAL_FRAME!r}; got {frame!r}')
