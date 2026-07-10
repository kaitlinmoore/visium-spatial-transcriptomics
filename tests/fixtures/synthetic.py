"""Deterministic synthetic Visium-like AnnData for the test suite.

This is the committed fixture, expressed as reproducible code rather than a
binary ``.h5ad`` blob: it is reviewable, diff-able, and regenerable, which suits
a learning project whose whole premise is explaining what the data is.

Geometry
--------
Spots sit on a real hexagonal lattice in *full-resolution pixel* coordinates
(``obsm["spatial"]``), the same frame a true Visium section uses. On a triangular
lattice every interior spot has exactly 6 equidistant neighbors — the Visium hex
neighborhood — so ``sq.gr.spatial_neighbors(coord_type="grid", n_neighs=6)`` and
the ``build_graph`` bridge have something honest to recover, and
``assert_hex_neighbors`` has something real to assert. Pixel layout for a
pointy-top lattice with spacing ``s``::

    x = col * s + (row % 2) * (s / 2)
    y = row * (s * sqrt(3) / 2)

Genes (raw integer counts in ``X``)
-----------------------------------
- ``GENE_HH``   : a planted High-High cluster — high counts in one contiguous
                  patch, low elsewhere. Strong positive global Moran's I; its HH
                  LISA spots should be exactly the planted patch.
- ``GENE_DISP`` : dispersed / spatially random counts. Global Moran's I ~ 0.
- ``GENE_GRAD`` : a smooth left-to-right gradient. Positive autocorrelation of a
                  different shape than the patch, to exercise the ranking.

Everything is seeded, so the fixture is identical on every machine.
"""

from __future__ import annotations

import numpy as np

# --- lattice parameters ------------------------------------------------------
N_ROWS = 7
N_COLS = 7
SPACING = 100.0  # full-resolution pixels between adjacent spots (Visium ~100um)
ROW_HEIGHT = SPACING * np.sqrt(3.0) / 2.0

GENE_NAMES = ["GENE_HH", "GENE_DISP", "GENE_GRAD"]

# Center of the planted High-High patch, in (row, col) lattice indices, and the
# lattice radius (in spots) of the patch. Recorded so tests can assert recovery.
HH_CENTER_RC = (3, 3)
HH_RADIUS = 1  # -> the center spot plus its 6 hex neighbors


def _lattice():
    """Return (rows, cols, xy) for the full triangular lattice.

    ``xy`` is the (n, 2) array of full-resolution pixel coordinates.
    """
    rows, cols, xy = [], [], []
    for r in range(N_ROWS):
        for c in range(N_COLS):
            x = c * SPACING + (r % 2) * (SPACING / 2.0)
            y = r * ROW_HEIGHT
            rows.append(r)
            cols.append(c)
            xy.append((x, y))
    return np.array(rows), np.array(cols), np.array(xy, dtype=float)


def _hh_patch_mask(rows: np.ndarray, cols: np.ndarray, xy: np.ndarray) -> np.ndarray:
    """Boolean mask of spots inside the planted High-High patch.

    Defined geometrically (a small radius in pixel space around the patch
    center) so the patch is the true hex neighborhood, independent of indexing.
    """
    cr, cc = HH_CENTER_RC
    center_x = cc * SPACING + (cr % 2) * (SPACING / 2.0)
    center_y = cr * ROW_HEIGHT
    dist = np.hypot(xy[:, 0] - center_x, xy[:, 1] - center_y)
    # radius just past one spacing catches the center + its 6 neighbors
    return dist <= SPACING * (HH_RADIUS + 0.1)


def make_synthetic_visium(*, seed: int = 0):
    """Build the synthetic Visium-like AnnData.

    Returns
    -------
    anndata.AnnData
        ``X`` = raw integer counts (n_spots x 3), ``obsm["spatial"]`` = full-res
        pixel coords, ``obs`` carries ``array_row``/``array_col`` and the boolean
        ``in_hh_patch``, and ``uns`` records the coordinate frame + scalefactors.
    """
    import anndata as ad
    import pandas as pd

    rng = np.random.default_rng(seed)
    rows, cols, xy = _lattice()
    n = xy.shape[0]

    hh_mask = _hh_patch_mask(rows, cols, xy)

    # GENE_HH: high counts inside the patch, low background. Clear HH cluster.
    gene_hh = np.where(
        hh_mask,
        rng.poisson(20.0, size=n),
        rng.poisson(1.0, size=n),
    ).astype(np.float32)

    # GENE_DISP: spatially random counts, no location dependence.
    gene_disp = rng.poisson(5.0, size=n).astype(np.float32)

    # GENE_GRAD: smooth gradient along x (autocorrelated, non-patch shape).
    x_norm = (xy[:, 0] - xy[:, 0].min()) / (np.ptp(xy[:, 0]) + 1e-9)
    gene_grad = rng.poisson(2.0 + 18.0 * x_norm).astype(np.float32)

    X = np.column_stack([gene_hh, gene_disp, gene_grad]).astype(np.float32)

    obs = pd.DataFrame(
        {
            "array_row": rows.astype(int),
            "array_col": cols.astype(int),
            "in_hh_patch": hh_mask,
        },
        index=[f"spot_{i}" for i in range(n)],
    )
    var = pd.DataFrame(index=GENE_NAMES)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = xy
    adata.uns["spatial_frame"] = "fullres_pixels"
    adata.uns["synthetic"] = {
        "seed": seed,
        "spacing": SPACING,
        "hh_center_rc": HH_CENTER_RC,
        "hh_radius": HH_RADIUS,
    }
    return adata


def hh_patch_indices(*, seed: int = 0) -> np.ndarray:
    """Integer positions of the planted High-High spots (for recovery tests)."""
    rows, cols, xy = _lattice()
    return np.flatnonzero(_hh_patch_mask(rows, cols, xy))


if __name__ == "__main__":  # optional: write an .h5ad for manual inspection
    import pathlib

    out = pathlib.Path(__file__).with_name("synthetic_visium.h5ad")
    make_synthetic_visium().write_h5ad(out)
    print(f"wrote {out} (gitignored — regenerate with `python {pathlib.Path(__file__).name}`)")
