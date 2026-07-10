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
    xy = np.asarray(spatial, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2:
        raise ValueError(f"spatial must be (n, 2); got shape {xy.shape}")
    if not tissue_hires_scalef > 0:
        raise ValueError(f"tissue_hires_scalef must be > 0; got {tissue_hires_scalef}")
    return xy * float(tissue_hires_scalef)


def plot_hotspots(adata, values, *, scalefactor: float, image=None, ax=None, title: str = ""):
    """Scatter per-spot ``values`` (e.g. LISA quadrants or Gi* z) on the H&E.

    Coordinates are scaled via :func:`to_hires_coords` before plotting, so the
    scatter lands in the SAME hires-pixel frame as ``image`` (the anti-
    misalignment guarantee). ``values`` may be numeric (a colorbar is drawn) or
    categorical labels like LISA quadrants (a legend is drawn). ``image`` is the
    hires PNG array used purely as a background.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    values = np.asarray(values)
    if values.shape[0] != adata.n_obs:
        raise ValueError(
            f"values length {values.shape[0]} != n_obs {adata.n_obs}"
        )

    xy = to_hires_coords(adata.obsm["spatial"], scalefactor)
    if ax is None:
        _, ax = plt.subplots()

    # The image establishes the hires-pixel frame (origin upper, y downward); the
    # scaled spots live in that same frame, so no y-flip is applied on top of it.
    if image is not None:
        ax.imshow(image)

    if np.issubdtype(values.dtype, np.number):
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=values, cmap="coolwarm", s=20)
        ax.figure.colorbar(sc, ax=ax)
    else:
        from matplotlib.patches import Patch

        cats = list(dict.fromkeys(values.tolist()))
        cmap = plt.get_cmap("tab10")
        color_for = {c: cmap(i % 10) for i, c in enumerate(cats)}
        ax.scatter(xy[:, 0], xy[:, 1], c=[color_for[v] for v in values], s=20)
        ax.legend(handles=[Patch(color=color_for[c], label=str(c)) for c in cats])

    # Without a background image, invert y so spots keep the image pixel
    # convention (y increases downward) rather than matplotlib's default y-up.
    if image is None:
        ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_title(title)
    return ax
