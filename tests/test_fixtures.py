"""Sanity checks on the synthetic fixture itself.

These run today (no pipeline code required) and prove two things the rest of the
suite leans on: the coordinates form a genuine hex lattice, and the planted
High-High gene really is spatially concentrated in the patch. If these break, a
fixture change silently moved the ground truth, so catch it here.
"""

from __future__ import annotations

import numpy as np


def test_shape_and_genes(adata):
    assert adata.n_obs == 49  # 7 x 7 lattice
    assert list(adata.var_names) == ["GENE_HH", "GENE_DISP", "GENE_GRAD"]
    assert adata.obsm["spatial"].shape == (adata.n_obs, 2)
    assert adata.uns["spatial_frame"] == "fullres_pixels"


def test_interior_spot_has_six_equidistant_neighbors(adata):
    """A triangular lattice: interior spots have exactly 6 nearest neighbors at
    the same distance. This is the hex neighborhood build_graph must recover."""
    xy = adata.obsm["spatial"]
    d = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)

    nn = d.min(axis=1)  # nearest-neighbor distance per spot
    interior = []
    for i in range(adata.n_obs):
        # count neighbors within a hair of this spot's nearest distance
        close = np.isclose(d[i], nn[i], rtol=1e-6)
        if close.sum() == 6:
            interior.append(i)
    # the 7x7 lattice has a solid interior block with full hex neighborhoods
    assert len(interior) >= 9


def test_hh_patch_is_the_center_hex(adata, hh_indices):
    """The planted patch is a center spot plus its 6 neighbors (7 spots)."""
    assert len(hh_indices) == 7
    assert adata.obs["in_hh_patch"].to_numpy().sum() == 7


def test_planted_gene_is_concentrated(adata, hh_indices):
    """GENE_HH is far higher inside the patch than outside; GENE_DISP is not."""
    hh = adata[:, "GENE_HH"].X.ravel()
    disp = adata[:, "GENE_DISP"].X.ravel()
    mask = np.zeros(adata.n_obs, dtype=bool)
    mask[hh_indices] = True

    assert hh[mask].mean() > 5 * hh[~mask].mean()
    # dispersed gene should show no such patch contrast
    assert disp[mask].mean() < 2 * disp[~mask].mean()


def test_scalefactors_present(scalefactors):
    assert scalefactors["tissue_hires_scalef"] == 0.15
    assert 0 < scalefactors["tissue_lowres_scalef"] < scalefactors["tissue_hires_scalef"]
