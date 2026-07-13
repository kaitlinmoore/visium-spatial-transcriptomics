# visium-spatial-statistics

Global-to-local spatial autocorrelation of a human lymph node 10x Visium section.

## What this is

A spatial-statistics analysis of a single human lymph node Visium slide. The headline is spatial autocorrelation taken from the global statistic (Moran's I) down to **local** indicators (local Moran's I / LISA and Getis-Ord Gi\*) computed on the same spatial weights graph, with the resulting hotspots overlaid on the aligned H&E and checked against known immune-compartment architecture. Neighborhood enrichment across clusters is a secondary, tissue-architecture read-out.

This is an analysis / inference project, not a predictive one. The point is methodology: constructing the spatial weights graph correctly, making the global→local move that a standard spatial-transcriptomics pipeline stops short of, running permutation inference with disciplined multiple testing, and validating hotspots against tissue structure rather than asserting them.

The spatial statistics here are the same family used in GIS: the Visium spot graph *is* a spatial weights matrix, Moran's I is the same statistic, and the LISA quadrant interpretation is identical, using spots instead of polygons, a gene instead of a socioeconomic variable. The local indicators are the part squidpy does not provide out of the box. They come from the PySAL toolkit (`esda`/`libpysal`).

**Biological framing.** A lymph node has crisp, well-separated compartments (B-cell follicles and germinal centers, the T-cell paracortex, sinuses). If the spatial statistics are working, the High-High LISA clusters for compartment markers should recover that architecture, showing follicle markers over the follicles, T-zone markers over the paracortex. That recovery is the validation, and should be checked against literature markers with sources, not eyeballed.

## Data

Single Visium section. This is a coordinate grid, not a sequence-alignment problem. Raw data is gitignored. The table below is the provenance record, filled in at download.

| Field | Value |
|---|---|
| Dataset ID | `V1_Human_Lymph_Node` |
| Assay | 10x Visium (spatial gene expression) |
| Source | `sq.datasets.visium("V1_Human_Lymph_Node")` → https://exampledata.scverse.org/squidpy/10x_genomics/V1_Human_Lymph_Node/ (scverse mirror of 10x Genomics' "Human Lymph Node" section) |
| Space Ranger version | `spaceranger-1.1.0` (from the matrix `.h5` `software_version`); chemistry Spatial 3' v1 |
| Reference transcriptome | `GRCh38` (the 10x reference the counts were generated against — provenance only; the analysis does not use it) |
| Access date | 2026-07-10 |
| Spots | 4,035 raw → 4,025 post-QC (defaults: `min_counts=500`, `min_genes=250`, `max_pct_mito=30`) |
| Genes | 36,601 raw → 19,812 post-QC (`min_spots=10`) |
| Images | `tissue_hires_image.png`, `tissue_lowres_image.png`, `scalefactors_json.json` |
| Access date | `<...>` |

The committed synthetic fixture (`tests/fixtures/synthetic.py`, `make_synthetic_visium`) is wired in as a **smoke-test path** so the pipeline runs end to end on a fresh clone with **no download**, before the lymph node section is fetched (see *How to run*). To acquire the real section, `load_visium.load_visium_dataset("V1_Human_Lymph_Node")` downloads via squidpy into `data/visium/…` (gitignored) and re-reads it through the owned `load_visium`.

## Environment & setup

Runs under WSL2 (Ubuntu). Dependencies are uv-managed (`pyproject.toml` + committed `uv.lock`), Python 3.11.

```bash
uv sync
```

Core stack: `squidpy`, `scanpy`, `anndata` (assay + graph + enrichment), `esda`, `libpysal` (local spatial statistics), `leidenalg`, `igraph` (clustering), plus `numpy`/`pandas`/`matplotlib`/`pytest`. All from PyPI wheels; no conda.

### Project layout

```
visium-spatial-statistics/
  data/                        # gitignored: Space Ranger output (matrix, spatial/, images), derived AnnData
  src/
    visium_spatial/            # importable package (installed editable by `uv sync`)
      __init__.py              # pipeline map; imports stay lazy so squidpy/esda load only when used
      load_visium.py           # sq.read.visium wrapper + coordinate/scalefactor extraction and checks
      qc.py                    # spot/gene QC and filtering
      preprocess.py            # normalize_total + log1p + HVG (scanpy)
      cluster.py               # neighbors + Leiden (compartment proxy)
      moran_scratch.py         # from-scratch global Moran's I (owned core; validated vs squidpy + esda)
      build_graph.py           # sq.gr.spatial_neighbors (hex) + squidpy->libpysal weights bridge (owned core)
      global_autocorr.py       # sq.gr.spatial_autocorr wrapper -> ranked spatially variable genes
      local_autocorr.py        # esda LISA (Moran_Local) + Getis-Ord Gi* (G_Local) on the shared graph
      multitest.py             # FDR within a gene; isolate/island handling (owned core)
      overlay.py               # hotspot maps overlaid on the aligned H&E (scalefactor-aware)
  notebooks/
    eda.ipynb                  # QC, coordinate + scalefactor sanity, clustering (+ smoke-test cell)
    autocorr.ipynb             # global -> local autocorrelation, hotspot overlays
    nhood.ipynb                # neighborhood enrichment, compartment adjacency
  tests/                       # pytest against committed synthetic fixtures
    conftest.py                # shared fixtures over the synthetic AnnData
    fixtures/                  # tiny synthetic Visium-like AnnData builder + scalefactors.json
  docs/methodology.md
  README.md
  CLAUDE.md
  pyproject.toml               # uv-managed; hatchling build backend
  uv.lock                      # committed
  LICENSE                      # MIT
```

## How to run

**1. Smoke test — no download.** Each notebook's first section runs the real package code against the committed synthetic fixture (`make_synthetic_visium`), to confirm the environment and the owned code work on a fresh clone. In `notebooks/eda.ipynb` this builds the graph, asserts the hex neighborhood, and ranks by global Moran's I via `visium_spatial.build_graph` and `visium_spatial.global_autocorr`.

**2. Real analysis — after downloading the lymph node section** (see *Data* and the action items in the planning package):

```bash
uv run jupyter lab
```

Then run the notebooks in order:

1. `notebooks/eda.ipynb` — load via `sq.read.visium`, QC, and the coordinate + scalefactor sanity checks (the pre-check that the overlay will align).
2. `notebooks/autocorr.ipynb` — global Moran's I ranking (validated against the from-scratch implementation), then LISA + Gi\* on the top-ranked genes, FDR correction, and hotspot overlays on the aligned H&E with the compartment-recovery check.
3. `notebooks/nhood.ipynb` — Leiden clustering as a compartment proxy, then neighborhood enrichment across clusters.

The `visium_spatial` package is installed editable by `uv sync`, so the notebooks just `from visium_spatial.build_graph import …`; the notebooks orchestrate the modules.

## Method summary

- **One shared weights matrix.** The Visium spatial graph is built once with `sq.gr.spatial_neighbors` on the hexagonal grid (`coord_type="grid"`, `n_neighs=6`) and bridged into a `libpysal` weights object, so the global (squidpy) and local (esda) statistics run on the *same* matrix. Row-standardized for analytic p-values.
- **Global as the cheap gate, local as the deliverable.** Genes are ranked by global Moran's I first; LISA/Gi\* run only on the top set. Cheap signal first, the expensive per-spot pass conditioned on hits.
- **Owned statistic, validated.** A from-scratch global Moran's I (weights matrix, row-standardization, the I formula, a permutation null) is the whiteboard-defensible core; it is checked against both squidpy and esda on the same weights before the local layer is trusted.
- **Multiple testing and isolates are first-class.** Per-spot p-values are FDR-corrected within a gene; comparison across genes is treated cautiously. Spots with no neighbors (tissue-boundary isolates) are handled explicitly and reported, never dropped silently.
- **Image as an aligned backdrop.** Hotspots are overlaid on the hires H&E by applying `tissue_hires_scalef` to the full-resolution pixel coordinates in `obsm["spatial"]`. No image-feature extraction — the histology is a registered background, not a data source.

Full rationale, the spatial-graph-as-weights-matrix crossover, the validation, and the coordinate/scalefactor hygiene are in [`docs/methodology.md`](docs/methodology.md).

## Results

Run on `V1_Human_Lymph_Node` (Space Ranger 1.1.0, GRCh38), seed 0. Reproduce via `notebooks/autocorr.ipynb`.

- **Data & QC.** 4,035 → 4,025 spots and 36,601 → 19,812 genes after QC (`min_counts=500`, `min_genes=250`, `max_pct_mito=30`, `min_spots=10`); 10 duplicate gene symbols disambiguated; 3 boundary spots are isolates (handled, not dropped).
- **Global validation.** The from-scratch Moran's I reproduces squidpy and esda on the shared row-standardized weights (e.g. `GENE_GRAD` synthetic parity `0.637977`); on the section, the Grasso-panel compartment markers rank as strongly spatially structured (`IGHG1` I=0.73, `FDCSP` 0.70, `CR2` 0.58, `TRBC1` 0.36; all q≈0).
- **Compartment recovery (the headline).** Using the published marker→domain panel (Grasso et al. 2025, *Eur. J. Immunol.*), High-High LISA hotspots for the four compact compartments form a clean block structure: **mean Jaccard 0.49 within vs 0.02 between** (follicle, germinal center, T-zone, medulla), each block completely segregated (0.00), with germinal centers nested in follicles — matching known architecture. Thin structures (blood/lymphatic vessels, the B–T interface) do not form hotspots — an honest limit of a compartment-hotspot method. See `docs/methodology.md` §8.
- **External benchmark (non-circular, the strongest result).** Scored against `germ_center`, a morphology-drawn germinal-centre annotation (368 spots, cell2location/Chrysalis processing of this section) that is spatially contiguous and *not* an abundance threshold: our GC-marker LISA hotspots recover it at **combined AUC 0.925**, and the BCL6/MYBL1 hotspot set is **92% precise** (154/168 spots inside annotated GCs, odds ratio 185) at 0.42 recall (FDR-gating captures GC cores). An expression-derived spatial statistic recovering an independently drawn morphological annotation.
- **Hotspot overlays** on the aligned H&E — representative markers (`FDCSP` follicle vs `TRBC1` T-zone) and our GC hotspots vs the annotated germinal centres — are produced in `autocorr.ipynb`.
- **Neighborhood enrichment (corroborating read-out).** 11 Leiden clusters; the B-like (`MS4A1`) and T-like (`CD3D`) clusters are each self-enriched (spatial-neighborhood z ≈ +29 / +70) and mutually depleted (cross z ≈ −11.6) — an independent method recovering the same follicle-vs-paracortex segregation. See `nhood.ipynb`.

## Limitations & scope

- **Single section.** Cross-section reproducibility (running a second lymph node slide and checking whether the same-marker hotspots reproduce) is a stretch, not part of v1.
- **Spot-level resolution.** Visium spots are multi-cell; nothing here is single-cell, and "compartment" is a spatial-domain notion, not a segmentation.
- **Recovery is validated, with caveats.** Beyond internal marker concordance, the GC hotspots are benchmarked against an external morphology-drawn annotation (AUC 0.925). But that annotation still correlates with GC expression (not a fully independent segmentation), recall is moderate (0.42 — GC cores, not full extents), and it is a single section.
- **No image analysis in v1.** The H&E is a backdrop; morphology features are a deferred stretch.

## License

MIT.
