# Methodology

> Working document. It records *why* each transform, statistic, and graph choice
> was made — the alternatives considered and why this one won — per the
> "explain as you build" principle. Sections are stubs to fill in as the
> pipeline lands; the reasoning already committed in the module docstrings is the
> source of truth until then.

## 1. Data and provenance
- Dataset, source, Space Ranger version, reference transcriptome, access date:
  see the provenance table in the README. Single section; no genome build, no
  CHROM:POS join.

## 2. Coordinate frame & scalefactor hygiene
- `obsm["spatial"]` is full-resolution image pixel coordinates (x, y). Which
  frame is in use is recorded, not assumed.
- Array (row/col) vs pixel (x/y) and y-axis orientation differ between image and
  plot conventions — state the convention wherever coordinates are drawn.
- Overlays multiply full-res coordinates by `tissue_hires_scalef`; mixing
  full-res coordinates with the scaled image is a silent misalignment (the
  tissue-image analog of a CRS/projection mismatch).

## 3. The spatial weights graph (the load-bearing choice)
- One graph, built once with `sq.gr.spatial_neighbors(coord_type="grid",
  n_neighs=6)` — the Visium hex neighborhood — and shared by the global
  (squidpy) and local (esda) statistics via a `libpysal` bridge.
- Neighbor structure is **asserted**, not assumed (interior degree = 6, isolates
  counted).
- Row-standardized (`transform="r"`) for analytic p-values.
- _Alternatives:_ `coord_type="generic"` (Delaunay / radius) — rejected, ignores
  known hex geometry and can bridge tissue gaps.

## 4. Global Moran's I — owned implementation and validation
- From-scratch statistic (weights, row-standardization, the I formula, a
  permutation null) is the defensibility core.
- It must reproduce squidpy's and esda's global Moran's I on the same weights,
  within tolerance, before the local layer is trusted.
- _Alternatives:_ analytic randomization z-scores (kept as cross-check only);
  Geary's C (measures local differences, not the global covariance).

## 5. Global → local
- Rank genes by global Moran's I first (cheap gate); run LISA / Gi* only on the
  top set (expensive per-spot pass).

## 6. Local indicators — LISA and Getis-Ord Gi*
- `esda.Moran_Local` (LISA) gives the HH/LH/LL/HL quadrant typology; High-High
  clusters for compartment markers are the compartment signal.
- `esda.G_Local(star=True)` (Gi*) gives a complementary hot/cold z-score.
- Both on the same shared weights.

## 7. Multiple testing & isolates
- Per-spot p-values FDR-corrected (Benjamini-Hochberg) **within a gene**;
  across-gene comparison treated cautiously and flagged.
- Isolates (no-neighbor boundary spots) handled explicitly (`island_weight`),
  counted, and reported — never dropped silently.
- _Alternatives:_ Bonferroni (too conservative); dropping isolates (dishonest —
  silently changes `n` and the tissue footprint).

## 8. Validation against tissue architecture
- Check whether High-High LISA clusters for known markers recover known immune
  compartments (B-cell follicle markers → follicles; T-zone markers →
  paracortex). Cite marker sources; qualitative validation, not a benchmark.

## 9. Secondary read-out: neighborhood enrichment
- Leiden clusters (expression graph — distinct from the spatial graph) as a
  compartment proxy, then `sq.gr.nhood_enrichment` across clusters.

## 10. Limitations
- Single section; spot-level (multi-cell) resolution; recovery is qualitative;
  no image-feature analysis in v1. See the README limitations section.
