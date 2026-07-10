"""Spot- and gene-level quality control and filtering.

Standard scanpy QC, kept separate from loading and from normalization so the
thresholds are visible and auditable rather than buried in a mega-function.

- Spot QC: total counts, number of genes per spot, mitochondrial fraction.
  Low-count / high-mito spots are tissue-edge or damaged and distort the
  spatial signal.
- Gene QC: drop genes detected in too few spots (they carry no spatial
  information and only inflate the multiple-testing burden downstream).

Choice: filter genes *before* ranking by global Moran's I, so the FDR
denominator reflects genes that could plausibly show structure. Alternative
(rank everything, filter after) inflates the correction and wastes the per-spot
local pass on noise.
"""

from __future__ import annotations


def compute_qc_metrics(adata):
    """Annotate ``adata`` with per-spot and per-gene QC metrics (in place).

    Wraps ``scanpy.pp.calculate_qc_metrics`` plus a mitochondrial flag/fraction
    (genes named ``MT-*``). ``percent_top=None`` so it works regardless of how
    many genes are present.
    """
    import scanpy as sc

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
    return adata


def filter_spots(adata, *, min_counts: int = 500, min_genes: int = 250, max_pct_mito: float = 30.0):
    """Return a copy of ``adata`` with low-quality spots removed.

    Thresholds are explicit arguments so they land in the notebook record, not a
    hidden default. Metrics are computed first if absent.
    """
    if "total_counts" not in adata.obs or "pct_counts_mt" not in adata.obs:
        compute_qc_metrics(adata)

    keep = (
        (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mito)
    )
    return adata[keep.values].copy()


def filter_genes(adata, *, min_spots: int = 10):
    """Drop genes detected in fewer than ``min_spots`` spots (in place)."""
    import scanpy as sc

    sc.pp.filter_genes(adata, min_cells=min_spots)
    return adata
