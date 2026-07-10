"""visium_spatial — spatial-statistics analysis of a 10x Visium section.

Pipeline (see docs/methodology.md and CLAUDE.md for the full rationale):

    load_visium      load the Visium AnnData + scalefactors, fix the coord frame
    qc               spot/gene QC and filtering
    preprocess       normalize_total + log1p + HVG
    cluster          expression-space neighbors + Leiden (compartment proxy)
    build_graph      build ONE hex spatial graph; bridge it into a libpysal W
                     shared by the global and local statistics
    moran_scratch    from-scratch global Moran's I — the owned, defensible check
    global_autocorr  rank spatially variable genes by global Moran's I (squidpy)
    local_autocorr   LISA (local Moran's I) + Getis-Ord Gi* on the top genes
    multitest        FDR within a gene; explicit isolate/island handling
    overlay          draw hotspots on the aligned H&E via the tissue scalefactor

Submodules are imported explicitly (e.g. ``from visium_spatial.build_graph
import build_spatial_graph``) so importing the package does not pull in the heavy
squidpy/esda stack until a specific stage is used.
"""

__version__ = "0.0.0"
