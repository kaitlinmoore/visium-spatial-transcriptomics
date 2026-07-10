"""Overlay hotspot maps on the aligned H&E, scalefactor-aware.

The histology is a **registered backdrop**, not a data source — no image-feature
extraction (that is a deferred stretch goal). This module's only job is to place
the spatial statistics on top of ``tissue_hires_image.png`` in the right spot.

The one hazard it exists to prevent: ``obsm["spatial"]`` is in full-resolution
pixel coordinates, but ``tissue_hires_image.png`` is downscaled. To draw a spot
on the hires image you must multiply its coordinates by ``tissue_hires_scalef``
from ``scalefactors_json.json``. Mixing full-res coordinates with the scaled
image is a silent misalignment — the tissue-image analog of a CRS/projection
mismatch — so the factor is applied explicitly here and never eyeballed.

Alternative considered: overlay on the lowres image. Same math with
``tissue_lowres_scalef``; hires is the default for legible hotspot maps.
"""

from __future__ import annotations

import numpy as np


def to_hires_coords(spatial: np.ndarray, tissue_hires_scalef: float) -> np.ndarray:
    """Scale full-resolution pixel coordinates onto the hires image.

    ``hires_xy = spatial_xy * tissue_hires_scalef``. The single explicit place
    the scalefactor is applied; everything else consumes the result.
    """
    raise NotImplementedError


def plot_hotspots(adata, values, *, scalefactor: float, image=None, ax=None, title: str = ""):
    """Scatter per-spot ``values`` (e.g. LISA quadrants or Gi* z) on the H&E.

    Coordinates are scaled via :func:`to_hires_coords` before plotting. ``image``
    is the hires PNG array used purely as a background.

    Returns
    -------
    matplotlib.axes.Axes
    """
    raise NotImplementedError
