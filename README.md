# Peri-implantitis mechanotransduction — reproducible analysis code

This repository contains the R and Python code used for the bioinformatic
analysis reported in the manuscript **"Integrated bulk and single-cell
transcriptomic analysis identifies mechanotransduction-related candidate
genes in peri-implantitis"**. The three candidate genes highlighted in the
study are **ICAM1, HPN, and COL19A1**.

The analysis is a secondary analysis of publicly deposited transcriptomic
data and does not include new experiments.

## Overview of the workflow

1. **Bulk discovery cohort (GSE223924):** TMM normalization and quality
   control, exploratory PCA, and limma-voom differential expression.
2. **Weighted gene co-expression network analysis (WGCNA):** module
   eigengenes were correlated with binary disease status; the disease-related
   red module was retained for candidate discovery.
3. **Three-way intersection:** WGCNA red-module genes were intersected with
   differentially expressed genes and a pre-existing curated
   mechanotransduction-related gene list (2,536 symbols).
4. **Enrichment, protein–protein interaction, and repeated-subsampling
   feature selection** (LASSO, SVM-RFE, random forest; 100 stratified
   subsamples) to prioritize candidate genes.
5. **Independent cohort assessment (GSE33774):** individual-gene expression,
   direction-adaptive ROC curves, and nominal Wilcoxon tests.
6. **Single-cell analysis (GSE272774):** Seurat processing, marker-based
   annotation, CellChat communication, virtual gene perturbation, trajectory
   analysis, and gene visualization.
7. **Complementary exploratory analyses:** immune deconvolution (CIBERSORT),
   GSEA/GSVA, consensus sample subtyping, compound–gene network, and
   molecular docking.

## Repository structure

```
.
├── config.R               # Auto-locates the repository root; defines paths
├── LICENSE                # MIT License
├── README.md
├── 17.venn/
│   └── MRG.txt            # Curated mechanotransduction gene list (2,536)
├── 27.GSEA/               # Place the MSigDB KEGG GMT here (see below)
├── 39.drugEnrich/         # Place the DSigDB reference table here (see below)
└── YM/                    # Analysis scripts, grouped by numbered step
    ├── 01.设计与简介/      # Package installation helpers
    ├── 03.bulk数据下载/    # Discovery cohort input check
    ├── 05.bulk标准化质控/  # Normalization and quality control
    ├── 07.箱线图/ ... 09.差异分析/
    ├── 10.单细胞数据下载/  # Single-cell download helpers (R/Python/PowerShell)
    ├── 12.Seurat分析/      # Seurat build, QC, integration, annotation
    ├── 13.细胞通讯/ ... 15.ClusterGVis/
    ├── 16.WGCNA/ ... 24.ROC/
    ├── 25.验证集差异/      # Independent cohort assessment
    ├── 27.GSEA/ ... 39.药物调控网络/
    └── 41.分子对接/        # Molecular docking instructions
```

Numbered directories follow the execution order. Scripts whose names begin
with `fix_` or `debug_` are reproducibility/repair scripts that document how
intermediate outputs were corrected; they are not part of the main pipeline.

## Requirements

The analysis was developed with **R 4.6.0** on Windows. A recent R (>= 4.3)
with the packages below should work. Key R packages:

- **Bioconductor:** edgeR, limma, clusterProfiler, org.Hs.eg.db,
  ComplexHeatmap, enrichplot, GSVA, GSEABase, slingshot, scTenifoldKnk,
  harmony
- **CRAN:** Seurat, WGCNA, ggplot2, ggrepel, pheatmap, cowplot, patchwork,
  dplyr, tidyr, tidyverse, Matrix, reticulate, igraph, httr, glmnet, e1071,
  randomForest, pROC, rms, ggpubr, corrplot, RColorBrewer, ggsci, circlize,
  ggvenn, ConsensusClusterPlus, ggplotify
- **GitHub:** CellChat (`jinworks/CellChat`), ClusterGVis, scRNAtoolVis and
  jjAnno

Use the installation helpers in `YM/01.设计与简介/` or install packages with
`BiocManager::install(...)` and `remotes::install_github(...)`.

The Python helper scripts use **Python 3** with `anndata`, `numpy`, `pandas`,
and `scanpy`. The R reticulate steps additionally require a Python environment
with `anndata`.

## Data setup

All source expression data are publicly available through GEO.

| Data | Accession | Place in repository |
|---|---|---|
| Discovery bulk RNA-seq counts | GSE223924 | `YM/03.bulk数据下载/` |
| Independent microarray cohort | GSE33774 | `YM/25.验证集差异/` (with GPL6244 annotation) |
| Single-cell gene–cell matrices | GSE272774 | `YM/10.单细胞数据下载/GSE272774/` |
| Integrated periodontal atlas (h5ad) | CELLxGENE collection | `YM/10.单细胞数据下载/periodontal_atlas.h5ad` |

CELLxGENE collection:
`https://cellxgene.cziscience.com/collections/71f4bccf-53d4-4c12-9e80-e73bfb89e398`
See `YM/10.单细胞数据下载/手动下载指引.md` for the manual download steps and
`download_zenodo.ps1` for a resumable Zenodo download helper.

Two third-party reference files are **not redistributed here** for licensing
and size reasons:

1. **MSigDB KEGG GMT** (`c2.cp.kegg.Hs.symbols.gmt`): download from the
   Molecular Signatures Database (free registration) and place it in
   `27.GSEA/`.
2. **DSigDB reference table** (`DSigDB_All_detailed.txt`, ~52 MB): download
   from the Drug Signatures Database and place it in `39.drugEnrich/`.

Large binaries and generated outputs (`.rds`, `.h5ad`, matrices, figures,
tables) are excluded via `.gitignore`; they are regenerated by running the
scripts in order. The final annotated Seurat object (`seu_annotated.rds`) is
rebuilt from GSE272774 and is not stored in this repository.

## Running the analysis

Each script automatically sources `config.R`, which locates the repository
root from the working directory, so no absolute paths need to be edited. Run
scripts from the repository (or from any subdirectory within it), for example:

```r
Rscript YM/05.bulk标准化质控/normalize.R
Rscript YM/09.差异分析/diff.R
Rscript YM/16.WGCNA/wgcna.R
```

Recommended order by numbered directory:

```
01 (install packages)
03 → 05 → 07 → 08 → 09          # bulk import, QC, PCA, differential expression
10 → 12                         # single-cell download and Seurat processing
13 → 14 → 15                    # CellChat and single-cell visualization
16 → 17 → 18 → 19               # WGCNA, intersection, enrichment, PPI
20 → 21 → 22 → 23 → 24          # feature selection, core genes, ROC, nomogram
25                              # independent cohort assessment
27 → 28 → 29 → 30 → 31          # GSEA, GSVA, CIBERSORT, immune analyses
32 → 33 → 34 → 35               # consensus subtyping and gene visualization
36 → 37 → 38 → 39               # trajectory, virtual knockout, CellChat, drugs
41                              # molecular docking (manual, see instructions)
```

Intermediate files written by an upstream step are read by the downstream
steps, so the order should be preserved on a first run.

## Methodological notes

- **Direction-adaptive ROC.** Individual-gene ROC curves use pROC with
  automatic direction selection; the predictor direction is chosen within each
  cohort. The reported AUCs therefore quantify direction-agnostic
  discrimination and do not test a prespecified biological direction or
  validate a fixed prediction rule.
- **Discovery AUCs.** The apparent AUC of 1.000 in the discovery cohort is an
  in-sample estimate from the same data used for feature selection and is
  expected to be optimistic.
- **Single-cell units.** Cells are not treated as independent donor replicates;
  pooled cell proportions and communication comparisons are descriptive.
- **Virtual perturbation.** scTenifoldKnk knock-out results are computational,
  network-based simulations with small effect sizes and are not experimental
  evidence.
- **Mechanotransduction framing.** The curated list is used as an annotation
  filter; this does not measure tissue loading or test cellular responses to
  mechanical force.

## Citation

Please cite the associated manuscript. Citation details and the full reference
list are provided in the manuscript. Source datasets: GSE223924, GSE33774, and
GSE272774.

## License

This code is released under the MIT License. See `LICENSE`; update the
copyright holder line with the author name before publication.
