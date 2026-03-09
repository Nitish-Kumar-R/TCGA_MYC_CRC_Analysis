# MYC-Centered Integrative Analysis of Colorectal Cancer (TCGA COAD/READ)

> End-to-end reproducible pipeline for MYC-driven colorectal cancer transcriptomics using TCGA COAD + READ RNA-seq data, with downstream machine learning survival prediction.

---

## Overview

This repository contains two complementary analysis pipelines:

1. **`TCGA_MYC_CRC_Pipeline.R`** — A comprehensive R pipeline that downloads and processes TCGA COAD/READ RNA-seq data (STAR counts), performs differential expression, WGCNA co-expression network analysis, GSEA, survival analysis, immune infiltration estimation, and builds a ceRNA regulatory network — all centered on MYC activity in colorectal cancer (CRC).

2. **`ML_XGBoost_SHAP.py`** — A Python machine learning module that trains an XGBoost classifier to predict 3-year mortality using WGCNA module gene expression, with SHAP-based model explainability.

---

## Repository Structure

```
.
├── TCGA_MYC_CRC_Pipeline.R     # Main R analysis pipeline (25 sections)
├── ML_XGBoost_SHAP.py          # XGBoost + SHAP ML module (ME4 module)
├── requirements.txt            # Python dependencies
├── figures/                    # Output: 19 publication-ready TIFF figures (600 DPI)
├── tables/                     # Output: 12 CSV result tables
└── docs/
    └── pipeline_overview.md    # Detailed section-by-section description
```

---

## Pipeline Sections (R)

| Section | Description | Output |
|---------|-------------|--------|
| 00 | Configuration & thread limits | — |
| 01 | Dependency installation (CRAN + Bioconductor) | — |
| 02 | Library loading | — |
| 03 | Helper functions | — |
| 04 | TCGA download & preparation | Raw counts |
| 05 | Sample annotation & QC | — |
| 06 | Gene filtering | — |
| 07 | DESeq2 – Tumor vs Normal | Table 01–02 |
| 08 | PCA | Figure 1 |
| 09 | Heatmap – Tumor vs Normal | Figure 2 |
| 10 | MYC stratification | — |
| 11 | DESeq2 – MYC-High vs MYC-Low | Figure 3, Table 03 |
| 12 | Heatmap – MYC-High vs MYC-Low | Figure 3 |
| 13 | Volcano plots | Figures 4A & 4B |
| 14 | GO enrichment – MYC-High upregulated | Figure 5, Table 04 |
| 15 | WGCNA – Tumor samples | Figures 6 & 7 |
| 16 | MYC-correlated module & hub genes | Figures 8, 9A & 9B, Tables 05–07 |
| 17 | GO enrichment – ME4 hub genes | Figure 10, Table 08 |
| 18 | GSEA – MYC + Hallmark gene sets | Figures 11A, 11B & 12, Tables 09–10 |
| 19 | Survival analysis – ME4 score KM + Cox | Figure 13 |
| 20 | ML dataset export (leakage-free split) | Figure 14, CSV for Python |
| 21 | Immune infiltration – ssGSEA | Figures 15 & 16 |
| 22 | Module vs immune correlation | Figure 17 |
| 23 | Immune checkpoint analysis | Figure 18 |
| 24 | MYC–Immune interaction network | Figure 19 |
| 25 | ceRNA network – miRTarBase + ENCORI | Tables 11 & 12 |

---

## Outputs

### Figures (19 total, 600 DPI TIFF)
- **Figure 1**: PCA – Tumor vs Normal
- **Figure 2**: Heatmap – Tumor vs Normal DE genes
- **Figure 3**: Heatmap – MYC-High vs MYC-Low
- **Figures 4A/4B**: Volcano plots
- **Figure 5**: GO enrichment – MYC-High upregulated genes
- **Figures 6–7**: WGCNA dendrogram & module–trait relationships
- **Figures 8, 9A/9B**: MYC-correlated module hub genes
- **Figure 10**: GO enrichment – ME4 module
- **Figures 11A/11B, 12**: GSEA enrichment plots
- **Figure 13**: Kaplan–Meier + Cox – ME4 survival score
- **Figure 14**: ML train/test split summary
- **Figures 15–16**: ssGSEA immune infiltration
- **Figure 17**: Module–immune cell correlation heatmap
- **Figure 18**: Immune checkpoint expression (MYC-High vs Low)
- **Figure 19**: MYC–Immune interaction network
- **Fig_ML_ROC_3yr_ME4_600dpi.tiff**: XGBoost ROC curve
- **Fig_ML_SHAP_3yr_ME4_600dpi.tiff**: SHAP beeswarm summary

### Tables (12 total, CSV)
- Tables 01–02: Tumor vs Normal DE results
- Table 03: MYC-High vs MYC-Low DE results
- Table 04: GO enrichment – MYC-High
- Tables 05–07: WGCNA hub genes & module membership
- Table 08: GO enrichment – ME4 module
- Tables 09–10: GSEA results
- Tables 11–12: ceRNA network edges & nodes

---

## Requirements

### R (≥ 4.5)

Dependencies are auto-installed by the pipeline (Section 01). Key packages:

**Bioconductor (v3.22):**
- `TCGAbiolinks`, `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`
- `WGCNA`, `GSVA`, `msigdbr`, `enrichplot`
- `ComplexHeatmap`, `pheatmap`, `circlize`

**CRAN:**
- `ggplot2`, `ggrepel`, `dplyr`, `survival`, `survminer`
- `igraph`, `ggraph`, `data.table`, `tibble`

### Python (≥ 3.9)

```bash
pip install -r requirements.txt
```

See `requirements.txt` for pinned versions.

---

## Usage

### Step 1 — Run the R Pipeline

```r
Rscript TCGA_MYC_CRC_Pipeline.R
```

This will:
- Auto-install all R/Bioconductor dependencies
- Download TCGA COAD + READ data via `TCGAbiolinks`
- Generate all figures in `figures/` and tables in `tables/`
- Export `ML_dataset_ME4_GeneSymbols_Survival.csv` for the Python module

> **Thread limit**: Defaults to 6 threads. Adjust `MAX_THREADS` in Section 00 as needed.

### Step 2 — Run the ML Module

```bash
python ML_XGBoost_SHAP.py
```

Requires `ML_dataset_ME4_GeneSymbols_Survival.csv` (exported by Section 17A of the R pipeline).

Outputs:
- `Fig_ML_ROC_3yr_ME4_600dpi.tiff`
- `Fig_ML_SHAP_3yr_ME4_600dpi.tiff`

---

## External Data (Section 25 only)

The ceRNA network analysis (Section 25) requires two manually downloaded files:

| File | Source |
|------|--------|
| `miRTarBase_MTI.csv` | https://mirtarbase.cuhk.edu.cn |
| `ENCORI_miRNA_lncRNA.tsv` | https://rnasysu.com/encori |

Place both files in the project root before running. Section 25 will be skipped gracefully if they are absent.

---

## Reproducibility

- R seed: `set.seed(123)` (Section 00)
- Python seed: `random_state=42` throughout
- Bioconductor version pinned to `3.19`
- Thread count capped at 6 across all parallel backends

---

## Citation

If you use this pipeline in your research, please cite the relevant tools:
- [TCGAbiolinks](https://doi.org/10.1093/nar/gkv1507)
- [DESeq2](https://doi.org/10.1186/s13059-014-0550-8)
- [WGCNA](https://doi.org/10.1186/1471-2105-9-559)
- [clusterProfiler](https://doi.org/10.1016/j.xinn.2021.100141)
- [GSVA](https://doi.org/10.1186/1471-2105-14-7)
- [XGBoost](https://doi.org/10.1145/2939672.2939785)
- [SHAP](https://proceedings.neurips.cc/paper/2017/hash/8a20a8621978632d76c43dfd28b67767-Abstract.html)

---

## Author

**Nitish Kumar R**

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
