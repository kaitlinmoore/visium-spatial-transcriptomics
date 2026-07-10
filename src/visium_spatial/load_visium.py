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
    raise NotImplementedError


def extract_scalefactors(adata) -> dict:
    """Return the Visium scalefactors as a plain ``dict`` of floats.

    Keys: ``spot_diameter_fullres``, ``tissue_hires_scalef``,
    ``tissue_lowres_scalef``, ``fiducial_diameter_fullres``.
    """
    raise NotImplementedError


def assert_spatial_frame(adata) -> None:
    """Fail loudly if ``obsm["spatial"]`` is missing or the frame is unrecorded.

    A guard, not a transform: catches the silent-bug class where coordinates are
    absent, in the wrong units, or of the wrong shape before any statistic runs.
    """
    raise NotImplementedError
