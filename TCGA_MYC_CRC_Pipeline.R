################################################################################
# TCGA COAD/READ – MYC-Centered Integrative Analysis
# Author  : Nitish Kumar R
#
# Description:
#   End-to-end reproducible pipeline for MYC-driven colorectal cancer
#   transcriptomics using TCGA COAD + READ RNA-seq data (STAR counts).
#
# Outputs:
#   figures/   – 19 publication-ready TIFF files (600 DPI)
#   tables/    – 12 CSV result tables
#
# Sections:
#   00. Configuration & thread limits
#   01. Dependency installation (CRAN + Bioconductor)
#   02. Library loading
#   03. Helper functions
#   04. TCGA download & preparation
#   05. Sample annotation & QC
#   06. Gene filtering
#   07. DESeq2 – Tumor vs Normal
#   08. PCA (Figure 1)
#   09. Heatmap – Tumor vs Normal (Figure 2)
#   10. MYC stratification
#   11. DESeq2 – MYC-High vs MYC-Low (Figure 3)
#   12. Heatmap – MYC-High vs MYC-Low (Figure 3)
#   13. Volcano plots (Figures 4A & 4B)
#   14. GO enrichment – MYC-High upregulated genes (Figure 5)
#   15. WGCNA – Tumor samples (Figures 6 & 7)
#   16. MYC-correlated module & hub genes (Figures 8, 9A & 9B)
#   17. GO enrichment – ME4 module hub genes (Figure 10)
#   18. GSEA – MYC + Hallmark gene sets (Figures 11A, 11B & 12)
#   19. Survival analysis – ME4 score KM + Cox (Figure 13)
#   20. ML dataset – leakage-free train/test split (Figure 14)
#   21. Immune infiltration – ssGSEA (Figures 15 & 16)
#   22. Module vs immune correlation (Figure 17)
#   23. Immune checkpoint analysis (Figure 18)
#   24. MYC–Immune interaction network (Figure 19)
#   25. ceRNA network – miRTarBase + ENCORI (Tables 11 & 12)
#
# External files needed only for Section 25:
#   miRTarBase_MTI.csv       https://mirtarbase.cuhk.edu.cn
#   ENCORI_miRNA_lncRNA.tsv  https://rnasysu.com/encori
################################################################################


################################################################################
## 00. CONFIGURATION & THREAD LIMITS
################################################################################

## ── Global seed & options ────────────────────────────────────────────────────
set.seed(123)
options(stringsAsFactors = FALSE)

## ── Thread / parallelism ceiling: 6 threads across ALL parallel backends ─────
MAX_THREADS <- 6L

# R-level parallelism (data.table, OpenBLAS, etc.)
options(Ncpus = MAX_THREADS)

# OpenMP thread cap (used by many compiled packages)
Sys.setenv(OMP_NUM_THREADS      = MAX_THREADS,
           OMP_THREAD_LIMIT     = MAX_THREADS,
           OPENBLAS_NUM_THREADS = MAX_THREADS,
           MKL_NUM_THREADS      = MAX_THREADS,
           VECLIB_MAXIMUM_THREADS = MAX_THREADS,
           NUMEXPR_NUM_THREADS  = MAX_THREADS)

# data.table thread cap
if (requireNamespace("data.table", quietly = TRUE))
  data.table::setDTthreads(MAX_THREADS)

## ── Output directories ───────────────────────────────────────────────────────
dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

message(sprintf("[00] Configuration complete – max threads: %d", MAX_THREADS))


################################################################################
## 01. DEPENDENCY INSTALLATION
################################################################################

## ── CRAN packages ────────────────────────────────────────────────────────────
cran_pkgs <- c(
  "ggplot2", "ggrepel", "dplyr", "stringr", "tidyr",
  "data.table", "tibble", "survival", "survminer", "scales",
  "igraph", "ggraph"
)

missing_cran <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(missing_cran) > 0) {
  message("[01] Installing missing CRAN packages: ",
          paste(missing_cran, collapse = ", "))
  install.packages(
    missing_cran,
    dependencies = TRUE,
    Ncpus        = MAX_THREADS,
    quiet        = TRUE
  )
}

## ── Bioconductor bootstrap ───────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

## ── Bioconductor packages ────────────────────────────────────────────────────
bioc_pkgs <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "apeglm",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "msigdbr",
  "enrichplot",
  "WGCNA",
  "pheatmap",
  "ComplexHeatmap",
  "circlize",
  "GSVA"
)

missing_bioc <- bioc_pkgs[!bioc_pkgs %in% rownames(installed.packages())]
if (length(missing_bioc) > 0) {
  message("[01] Installing missing Bioconductor packages: ",
          paste(missing_bioc, collapse = ", "))
  BiocManager::install(
    missing_bioc,
    ask    = FALSE,
    update = FALSE,
    Ncpus  = MAX_THREADS
  )
}

message("[01] All dependencies satisfied.")


################################################################################
## 02. LIBRARY LOADING
################################################################################

suppressPackageStartupMessages({
  # Bioconductor
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  library(enrichplot)
  library(WGCNA)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(GSVA)

  # CRAN
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(data.table)
  library(tibble)
  library(survival)
  library(survminer)
  library(scales)
  library(igraph)
  library(ggraph)
})

## ── Apply data.table thread cap (post-load) ──────────────────────────────────
data.table::setDTthreads(MAX_THREADS)

## ── WGCNA thread cap ─────────────────────────────────────────────────────────
## enableWGCNAThreads() uses multi-threading internally; cap to MAX_THREADS.
## This must be called AFTER library(WGCNA).
enableWGCNAThreads(nThreads = MAX_THREADS)

## ── CRITICAL: WGCNA vs stats::cor namespace fix ──────────────────────────────
## blockwiseModules() internally calls WGCNA::cor(), which accepts extra args
## (weights.x, weights.y, cosine) that stats::cor() does not.
## Convention used throughout this script:
##   cor <- WGCNA::cor   immediately BEFORE each blockwiseModules() call
##   cor <- stats::cor   immediately AFTER  each blockwiseModules() call
cor <- stats::cor   # safe default

message("[02] Libraries loaded.")


################################################################################
## 03. HELPER FUNCTIONS
################################################################################

## ── save_tiff_km: Kaplan–Meier plot to TIFF ──────────────────────────────────
save_tiff_km <- function(fit, data, filename, title,
                          legend_title, legend_labs, palette,
                          break_by = 1000, width = 7, height = 7) {
  tiff(filename, width = width, height = height,
       units = "in", res = 600, compression = "lzw")
  print(
    ggsurvplot(
      fit,
      data              = data,
      risk.table        = TRUE,
      pval              = TRUE,
      conf.int          = FALSE,
      censor            = TRUE,
      censor.size       = 2.5,
      censor.shape      = 124,
      risk.table.height = 0.28,
      break.time.by     = break_by,
      xlab              = "Time (days)",
      ylab              = "Overall Survival Probability",
      legend.title      = legend_title,
      legend.labs       = legend_labs,
      palette           = palette,
      ggtheme           = theme_classic(base_size = 14),
      title             = title
    )
  )
  dev.off()
  message("Saved: ", filename)
}

## ── ensembl_to_symbols: strip version suffix & map to HGNC symbols ───────────
ensembl_to_symbols <- function(expr_matrix) {
  ens_clean <- sub("\\..*$", "", rownames(expr_matrix))
  syms <- mapIds(
    org.Hs.eg.db,
    keys      = ens_clean,
    keytype   = "ENSEMBL",
    column    = "SYMBOL",
    multiVals = "first"
  )
  rownames(expr_matrix) <- syms
  expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), , drop = FALSE]
  expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), , drop = FALSE]
  expr_matrix
}

## ── make_volcano: publication-quality volcano plot ───────────────────────────
## Features:
##   • Full uncapped y-axis range
##   • Plain circle points (shape = 16) throughout
##   • Top 5 up + top 5 down gene labels via ggrepel
##   • Dashed reference lines at LFC = ±1 and -log10(0.05)
##   • theme_bw, legend on top, no minor gridlines
make_volcano <- function(df, lfc_col, padj_col, label_col,
                          color_col, color_pal,
                          up_level, down_level, title) {

  df <- df %>%
    mutate(
      padj_safe = ifelse(.data[[padj_col]] == 0,
                         .Machine$double.xmin,
                         .data[[padj_col]]),
      y_val = -log10(padj_safe)
    )

  top_up <- df %>%
    filter(.data[[color_col]] == up_level) %>%
    arrange(desc(abs(.data[[lfc_col]]))) %>%
    slice_head(n = 5)

  top_dn <- df %>%
    filter(.data[[color_col]] == down_level) %>%
    arrange(.data[[lfc_col]]) %>%
    slice_head(n = 5)

  labels <- bind_rows(top_up, top_dn)

  ggplot(df, aes(x = .data[[lfc_col]], y = y_val)) +
    geom_point(
      aes(color = .data[[color_col]]),
      shape = 16, size = 1.5, alpha = 0.65
    ) +
    geom_text_repel(
      data               = labels,
      aes(label          = .data[[label_col]]),
      size               = 3.2,
      box.padding        = 0.4,
      point.padding      = 0.3,
      segment.color      = "grey50",
      segment.size       = 0.3,
      min.segment.length = 0.2,
      max.overlaps       = Inf,
      force              = 4,
      direction          = "both",
      show.legend        = FALSE
    ) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", linewidth = 0.45, color = "grey35") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", linewidth = 0.45, color = "grey35") +
    scale_color_manual(values = color_pal) +
    scale_x_continuous(expand = expansion(mult = 0.05)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(
      title = title,
      x     = expression(Log[2] ~ "Fold Change"),
      y     = expression(-Log[10] ~ "Adjusted" ~ italic(p) * "-value"),
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin      = margin(12, 20, 10, 12),
      legend.position  = "top",
      legend.key.size  = unit(0.45, "cm"),
      legend.text      = element_text(size = 10),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 10, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3)
    )
}

## ── save_go_dotplot: GO enrichment dotplot to TIFF ───────────────────────────
save_go_dotplot <- function(ego_obj, filename, title,
                             show_n = 10, wrap_width = 38,
                             width = 8, height = 6) {
  p <- dotplot(ego_obj, showCategory = show_n, font.size = 11) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = wrap_width)) +
    scale_color_gradient(low  = "#2166AC", high = "#B2182B",
                         name = "Adjusted\np-value") +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_bw() +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin        = margin(20, 30, 15, 20),
      panel.grid.major.y = element_blank(),
      legend.position    = "right"
    ) +
    coord_cartesian(clip = "off")
  ggsave(filename, p, dpi = 600, width = width, height = height,
         device = "tiff", compression = "lzw")
  message("Saved: ", filename)
}

message("[03] Helper functions defined.")

################################################################################
## 04. TCGA DOWNLOAD & PREPARATION
################################################################################

dir.create("data", showWarnings = FALSE)

data_file <- "data/tcga_coad_read_rnaseq.rds"

if (file.exists(data_file)) {
  
  message("[04] Local TCGA dataset found – loading cached data …")
  
  rna_se <- readRDS(data_file)
  
} else {
  
  message("[04] Querying TCGA COAD + READ RNA-seq data …")
  
  query <- GDCquery(
    project       = c("TCGA-COAD", "TCGA-READ"),
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(
    query,
    method = "api",
    files.per.chunk = 20
  )
  
  rna_se <- GDCprepare(query)
  
  saveRDS(rna_se, data_file)
  
  message("[04] TCGA data downloaded and cached.")
}

counts  <- assay(rna_se)
coldata <- as.data.frame(colData(rna_se))

message(sprintf("[04] Raw matrix: %d genes × %d samples",
                nrow(counts), ncol(counts)))

################################################################################
## 05. SAMPLE ANNOTATION & QC
################################################################################

message("[05] Annotating samples …")
print(table(coldata$shortLetterCode))

coldata$condition <- dplyr::case_when(
  coldata$shortLetterCode == "TP" ~ "Tumor",
  coldata$shortLetterCode == "NT" ~ "Normal",
  TRUE                            ~ NA_character_
)

keep    <- !is.na(coldata$condition)
counts  <- counts[, keep]
coldata <- coldata[keep, ]
coldata$condition <- factor(coldata$condition, levels = c("Normal", "Tumor"))

message(sprintf("[05] Samples retained – Tumor: %d  Normal: %d",
                sum(coldata$condition == "Tumor"),
                sum(coldata$condition == "Normal")))


################################################################################
## 06. GENE FILTERING  (≥10 counts in ≥20 % of samples)
################################################################################

message("[06] Filtering low-count genes …")

keep_genes  <- rowSums(counts >= 10) >= 0.2 * ncol(counts)
counts_filt <- counts[keep_genes, ]

message(sprintf("[06] Genes retained: %d  (removed: %d)",
                nrow(counts_filt), sum(!keep_genes)))


################################################################################
## 07. DESEQ2 – TUMOR vs NORMAL
################################################################################

message("[07] Running DESeq2: Tumor vs Normal …")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = ~ condition
)
dds    <- DESeq(dds)
res_TN <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# Variance-stabilised expression (used in all downstream figures)
vsd       <- vst(dds, blind = TRUE)
norm_expr <- assay(vsd)

res_TN_df <- as.data.frame(res_TN) %>%
  rownames_to_column("ensembl") %>%
  mutate(
    gene_symbol = rowData(rna_se)$gene_name[
      match(ensembl, rownames(rowData(rna_se)))
    ],
    Significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE                              ~ "Not Significant"
    )
  )

write.csv(res_TN_df, "tables/01_DGE_Tumor_vs_Normal.csv", row.names = FALSE)

message(sprintf("[07] Tumor vs Normal – up: %d  down: %d",
                sum(res_TN_df$Significance == "Upregulated",   na.rm = TRUE),
                sum(res_TN_df$Significance == "Downregulated", na.rm = TRUE)))


################################################################################
## 08. FIGURE 1 – PCA: Tumor vs Normal
################################################################################

message("[08] Generating Figure 1 (PCA) …")

pca_df     <- plotPCA(vsd, intgroup = "condition", ntop = 500, returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c(Normal = "#00BFC4", Tumor = "#F8766D")) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggtitle("PCA – TCGA Colorectal Cancer: Tumor vs Normal") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/Figure_1_PCA_Tumor_vs_Normal.tiff",
       p_pca, dpi = 600, width = 6, height = 5,
       device = "tiff", compression = "lzw")
message("Saved: figures/Figure_1_PCA_Tumor_vs_Normal.tiff")


################################################################################
## 09. FIGURE 2 – HEATMAP: Top 50 DEGs, Tumor vs Normal
################################################################################

message("[09] Generating Figure 2 (Heatmap – Tumor vs Normal) …")

top50_ens <- res_TN_df %>%
  filter(!is.na(padj), !is.na(gene_symbol)) %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(ensembl)

expr_heat_TN <- norm_expr[top50_ens, ]
rownames(expr_heat_TN) <- res_TN_df$gene_symbol[
  match(top50_ens, res_TN_df$ensembl)
]

z_TN <- t(scale(t(expr_heat_TN)))
z_TN[z_TN >  2] <-  2
z_TN[z_TN < -2] <- -2

sample_order_TN <- rownames(coldata)[order(coldata$condition)]
ann_col_TN      <- data.frame(
  Condition  = coldata$condition,
  row.names  = rownames(coldata)
)
ann_col_TN      <- ann_col_TN[sample_order_TN, , drop = FALSE]
ann_colors_TN   <- list(Condition = c(Normal = "#00BFC4", Tumor = "#F8766D"))

tiff("figures/Figure_2_Heatmap_Tumor_vs_Normal.tiff",
     width = 8, height = 9, units = "in", res = 600, compression = "lzw")
pheatmap(
  z_TN[, sample_order_TN],
  annotation_col    = ann_col_TN,
  annotation_colors = ann_colors_TN,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 7,
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  breaks            = seq(-2, 2, length.out = 101),
  border_color      = NA,
  main              = "Top 50 DEGs: Tumor vs Normal (CRC)"
)
dev.off()
message("Saved: figures/Figure_2_Heatmap_Tumor_vs_Normal.tiff")


################################################################################
## 10. MYC STRATIFICATION
################################################################################

message("[10] Stratifying samples by MYC expression …")

myc_row <- which(rowData(rna_se)$gene_name == "MYC")
stopifnot("MYC not found in rowData$gene_name" = length(myc_row) == 1)

myc_exp <- norm_expr[myc_row, ]

coldata$MYC_group <- factor(
  ifelse(myc_exp >= median(myc_exp, na.rm = TRUE),
         "MYC_high", "MYC_low"),
  levels = c("MYC_low", "MYC_high")
)

message(sprintf("[10] MYC-high: %d  MYC-low: %d",
                sum(coldata$MYC_group == "MYC_high"),
                sum(coldata$MYC_group == "MYC_low")))


################################################################################
## 11. DESEQ2 – MYC-HIGH vs MYC-LOW  (with apeglm LFC shrinkage)
################################################################################

message("[11] Running DESeq2: MYC-High vs MYC-Low …")

dds_myc <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = ~ MYC_group
)
dds_myc$MYC_group <- factor(dds_myc$MYC_group,
                              levels = c("MYC_low", "MYC_high"))
dds_myc <- DESeq(dds_myc)

res_myc <- results(dds_myc,
                    contrast = c("MYC_group", "MYC_high", "MYC_low"))

res_myc <- lfcShrink(
  dds_myc,
  coef = "MYC_group_MYC_high_vs_MYC_low",
  res  = res_myc
)

res_myc_df <- as.data.frame(res_myc) %>%
  rownames_to_column("ensembl") %>%
  mutate(
    gene_symbol = rowData(rna_se)$gene_name[
      match(ensembl, rownames(rowData(rna_se)))
    ],
    Significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "MYC-High",
      padj < 0.05 & log2FoldChange < -1 ~ "MYC-Low",
      TRUE                              ~ "Not Significant"
    )
  ) %>%
  filter(!is.na(gene_symbol))

write.csv(res_myc_df, "tables/02_DGE_MYC_High_vs_Low.csv", row.names = FALSE)

message(sprintf("[11] MYC-High vs MYC-Low – up: %d  down: %d",
                sum(res_myc_df$Significance == "MYC-High", na.rm = TRUE),
                sum(res_myc_df$Significance == "MYC-Low",  na.rm = TRUE)))


################################################################################
## 12. FIGURE 3 – HEATMAP: Top 30 DEGs, MYC-High vs MYC-Low
################################################################################

message("[12] Generating Figure 3 (Heatmap – MYC-High vs MYC-Low) …")

tumor_mask <- coldata$condition == "Tumor"

top30_myc_ens <- res_myc_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 30) %>%
  pull(ensembl)

expr_heat_myc <- norm_expr[top30_myc_ens, tumor_mask]
rownames(expr_heat_myc) <- res_myc_df$gene_symbol[
  match(top30_myc_ens, res_myc_df$ensembl)
]

z_myc <- t(scale(t(expr_heat_myc)))
z_myc[z_myc >  2] <-  2
z_myc[z_myc < -2] <- -2

myc_sample_order <- rownames(coldata)[tumor_mask][
  order(coldata$MYC_group[tumor_mask])
]
ann_col_myc <- data.frame(
  MYC_Status = coldata$MYC_group[tumor_mask],
  row.names  = rownames(coldata)[tumor_mask]
)[myc_sample_order, , drop = FALSE]

ann_colors_myc <- list(
  MYC_Status = c(MYC_low = "#4575B4", MYC_high = "#D73027")
)

tiff("figures/Figure_3_Heatmap_MYC_High_vs_Low.tiff",
     width = 8, height = 8, units = "in", res = 600, compression = "lzw")
pheatmap(
  z_myc[, myc_sample_order],
  annotation_col    = ann_col_myc,
  annotation_colors = ann_colors_myc,
  show_colnames     = FALSE,
  show_rownames     = TRUE,
  fontsize_row      = 8,
  cluster_cols      = FALSE,
  cluster_rows      = TRUE,
  color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks            = seq(-2, 2, length.out = 101),
  border_color      = NA,
  main              = "Top 30 DEGs: MYC-High vs MYC-Low (CRC Tumors)"
)
dev.off()
message("Saved: figures/Figure_3_Heatmap_MYC_High_vs_Low.tiff")


################################################################################
## 13. FIGURES 4A & 4B – VOLCANO PLOTS
################################################################################

message("[13] Generating Figures 4A & 4B (Volcano plots) …")

## ── Figure 4A: Tumor vs Normal ───────────────────────────────────────────────
p_vol_TN <- make_volcano(
  df         = res_TN_df %>% filter(!is.na(padj)),
  lfc_col    = "log2FoldChange",
  padj_col   = "padj",
  label_col  = "gene_symbol",
  color_col  = "Significance",
  color_pal  = c(
    Upregulated       = "#D73027",
    Downregulated     = "#4575B4",
    `Not Significant` = "grey70"
  ),
  up_level   = "Upregulated",
  down_level = "Downregulated",
  title      = "Differential Expression: Tumor vs Normal (CRC)"
)
ggsave("figures/Figure_4A_Volcano_Tumor_vs_Normal.tiff",
       p_vol_TN, width = 7.5, height = 6, dpi = 600,
       device = "tiff", compression = "lzw")
message("Saved: figures/Figure_4A_Volcano_Tumor_vs_Normal.tiff")

## ── Figure 4B: MYC-High vs MYC-Low ──────────────────────────────────────────
p_vol_myc <- make_volcano(
  df         = res_myc_df %>% filter(!is.na(padj)),
  lfc_col    = "log2FoldChange",
  padj_col   = "padj",
  label_col  = "gene_symbol",
  color_col  = "Significance",
  color_pal  = c(
    `MYC-High`        = "#D73027",
    `MYC-Low`         = "#4575B4",
    `Not Significant` = "grey70"
  ),
  up_level   = "MYC-High",
  down_level = "MYC-Low",
  title      = "Differential Expression: MYC-High vs MYC-Low (CRC)"
)
ggsave("figures/Figure_4B_Volcano_MYC_High_vs_Low.tiff",
       p_vol_myc, width = 7.5, height = 6, dpi = 600,
       device = "tiff", compression = "lzw")
message("Saved: figures/Figure_4B_Volcano_MYC_High_vs_Low.tiff")


################################################################################
## 14. FIGURE 5 – GO ENRICHMENT: MYC-High Upregulated Genes
################################################################################

message("[14] Running GO enrichment for MYC-High upregulated genes …")

myc_up_genes <- res_myc_df %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  pull(gene_symbol)

ego_MYC <- enrichGO(
  gene          = myc_up_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_MYC_simp <- clusterProfiler::simplify(ego_MYC, cutoff = 0.7,
                                           by = "p.adjust", select_fun = min)

save_go_dotplot(
  ego_MYC_simp,
  filename = "figures/Figure_5_GO_BP_MYC_High_Upregulated.tiff",
  title    = "GO Biological Processes – MYC-High CRC"
)

write.csv(as.data.frame(ego_MYC_simp),
          "tables/03_GO_MYC_High_simplified.csv",
          row.names = FALSE)


################################################################################
## 15. FIGURES 6 & 7 – WGCNA: Soft-Threshold & Gene Dendrogram
#
#  Thread usage is capped at MAX_THREADS (6) via:
#    • enableWGCNAThreads(nThreads = MAX_THREADS) called in Section 02
#    • nThreads argument passed explicitly to blockwiseModules()
#    • cor namespace swap guards each blockwiseModules() call
################################################################################

message("[15] Running WGCNA on tumor samples …")

expr_tumor <- norm_expr[, tumor_mask]

## Select top 5,000 most variable genes by MAD
mad_vals <- apply(expr_tumor, 1, mad)
top5k    <- names(sort(mad_vals, decreasing = TRUE))[1:5000]
datExpr  <- t(expr_tumor[top5k, ])   # samples × genes

## WGCNA sample & gene QC
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message(sprintf("[15] WGCNA QC removed %d samples and %d genes",
                  sum(!gsg$goodSamples), sum(!gsg$goodGenes)))
}

## ── Figure 6: Soft-threshold selection ───────────────────────────────────────
powers <- 1:20
sft    <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)

tiff("figures/Figure_6_WGCNA_Soft_Threshold_Selection.tiff",
     width = 8, height = 4, units = "in", res = 600, compression = "lzw")
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = expression(R^2 ~ "Scale-Free Fit"),
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "blue", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers, col = "red")
dev.off()
message("Saved: figures/Figure_6_WGCNA_Soft_Threshold_Selection.tiff")

## ── blockwiseModules with 6-thread cap ───────────────────────────────────────
softPower <- 6

cor <- WGCNA::cor       # namespace swap IN
net <- blockwiseModules(
  datExpr,
  power             = softPower,
  nThreads          = MAX_THREADS,   # explicit 6-thread cap
  TOMType           = "unsigned",
  minModuleSize     = 50,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 0
)
cor <- stats::cor       # namespace swap OUT

moduleColors <- labels2colors(net$colors)
message(sprintf("[15] WGCNA modules identified: %d",
                length(unique(moduleColors))))

## ── Figure 7: Gene dendrogram & module colours ───────────────────────────────
tiff("figures/Figure_7_WGCNA_Gene_Dendrogram_Module_Colours.tiff",
     width = 10, height = 5, units = "in", res = 600, compression = "lzw")
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colours",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()
message("Saved: figures/Figure_7_WGCNA_Gene_Dendrogram_Module_Colours.tiff")


################################################################################
## 16. FIGURES 8, 9A & 9B – MYC-Correlated Module & Hub Genes
################################################################################

message("[16] Identifying MYC-correlated WGCNA module …")

myc_expr_tumor <- myc_exp[tumor_mask][rownames(datExpr)]
MEs            <- moduleEigengenes(datExpr, moduleColors)$eigengenes
moduleTraitCor <- WGCNA::cor(MEs, myc_expr_tumor, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

## ── Figure 8: Module–MYC correlation heatmap ─────────────────────────────────
tiff("figures/Figure_8_WGCNA_Module_MYC_Correlation_Heatmap.tiff",
     width = 4, height = 7, units = "in", res = 600, compression = "lzw")
labeledHeatmap(
  Matrix      = matrix(moduleTraitCor, ncol = 1),
  xLabels     = "MYC expression",
  yLabels     = rownames(moduleTraitCor),
  ySymbols    = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = signif(matrix(moduleTraitCor, ncol = 1), 2),
  main        = "Module–MYC Correlation"
)
dev.off()
message("Saved: figures/Figure_8_WGCNA_Module_MYC_Correlation_Heatmap.tiff")

## ── Identify top MYC-correlated module ───────────────────────────────────────
module_of_interest <- rownames(moduleTraitCor)[which.max(moduleTraitCor[, 1])]
target_color       <- gsub("^ME", "", module_of_interest)
ME4_genes          <- colnames(datExpr)[moduleColors == target_color]

message(sprintf("[16] MYC-correlated module: %s  colour: %s  genes: %d",
                module_of_interest, target_color, length(ME4_genes)))

## ── Hub gene table ────────────────────────────────────────────────────────────
kME              <- signedKME(datExpr, MEs)
geneSignificance <- cor(datExpr, myc_expr_tumor, use = "p")

ME4_col <- grep(paste0("kME", target_color), colnames(kME), value = TRUE)
stopifnot("ME4 kME column not found" = length(ME4_col) == 1)

hub_ME4 <- data.frame(
  Gene   = ME4_genes,
  kME    = kME[ME4_genes, ME4_col],
  GS_MYC = geneSignificance[ME4_genes, 1]
)

stopifnot(
  nrow(hub_ME4) == length(ME4_genes),
  !any(is.na(hub_ME4$kME)),
  !any(is.na(hub_ME4$GS_MYC))
)

hub_ME4$ensembl_clean <- sub("\\..*$", "", hub_ME4$Gene)
hub_ME4$GeneSymbol    <- mapIds(
  org.Hs.eg.db,
  keys      = hub_ME4$ensembl_clean,
  keytype   = "ENSEMBL",
  column    = "SYMBOL",
  multiVals = "first"
)
hub_ME4 <- hub_ME4 %>%
  filter(!is.na(GeneSymbol)) %>%
  arrange(desc(abs(kME)))

message(sprintf("[16] Annotated hub genes: %d", nrow(hub_ME4)))
message("Top 15 hub genes:")
print(head(hub_ME4[, c("GeneSymbol", "kME", "GS_MYC")], 15))

write.csv(data.frame(ensembl = ME4_genes),
          "tables/04_ME4_gene_list.csv", row.names = FALSE)
write.csv(hub_ME4,
          "tables/05_ME4_hub_genes_annotated.csv", row.names = FALSE)

ME4_eigengene        <- MEs[, module_of_interest]
names(ME4_eigengene) <- rownames(datExpr)

## ── Figure 9A: Top 50 MYC-correlated hub genes heatmap ───────────────────────
message("[16] Generating Figure 9A …")

top_hub_genes <- hub_ME4 %>%
  arrange(desc(abs(GS_MYC))) %>%
  slice_head(n = 50) %>%
  pull(Gene)

expr_top   <- norm_expr[top_hub_genes, tumor_mask]
rownames(expr_top) <- hub_ME4$GeneSymbol[match(top_hub_genes, hub_ME4$Gene)]

expr_top_z <- t(scale(t(expr_top)))
expr_top_z[expr_top_z >  4] <-  4
expr_top_z[expr_top_z < -4] <- -4

myc_annotation <- myc_exp[tumor_mask]
top_anno <- HeatmapAnnotation(
  MYC = myc_annotation,
  col = list(
    MYC = colorRamp2(
      c(min(myc_annotation), median(myc_annotation), max(myc_annotation)),
      c("#ccece6", "#66c2a4", "#238b45")
    )
  ),
  annotation_height   = unit(6, "mm"),
  annotation_name_gp  = gpar(fontsize = 9, fontface = "bold")
)

ht <- Heatmap(
  expr_top_z,
  name              = "Z-score",
  col               = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
  cluster_rows      = TRUE,
  cluster_columns   = TRUE,
  show_column_names = FALSE,
  show_row_names    = TRUE,
  row_names_gp      = gpar(fontsize = 8),
  column_title      = "Top 50 MYC-Correlated Hub Genes (WGCNA, CRC Tumors)",
  column_title_gp   = gpar(fontsize = 12, fontface = "bold"),
  top_annotation    = top_anno,
  heatmap_legend_param = list(
    title     = "Expression\n(Z-score)",
    at        = c(-4, -2, 0, 2, 4),
    title_gp  = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

tiff("figures/Figure_9A_WGCNA_Top_MYC_Correlated_Genes_Heatmap.tiff",
     width = 10, height = 9, units = "in", res = 600, compression = "lzw")
draw(ht)
dev.off()
message("Saved: figures/Figure_9A_WGCNA_Top_MYC_Correlated_Genes_Heatmap.tiff")

## ── Figure 9B: KEGG pathway enrichment of hub genes ──────────────────────────
message("[16] Running KEGG enrichment for hub genes (Figure 9B) …")

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys      = hub_ME4$GeneSymbol,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)

kegg_res <- enrichKEGG(
  gene         = entrez_ids,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

if (nrow(as.data.frame(kegg_res)) > 0) {
  kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID")
  p_kegg <- dotplot(kegg_res, showCategory = 12, font.size = 12) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    scale_color_gradient(low  = "#2166AC", high = "#B2182B",
                         name = "Adjusted\np-value") +
    labs(
      title = "KEGG Pathway Enrichment – MYC-Correlated Hub Genes (CRC)",
      x = "Gene Ratio", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.y = element_blank()
    )
  ggsave("figures/Figure_9B_KEGG_MYC_Correlated_Hub_Genes.tiff",
         p_kegg, width = 8, height = 6, dpi = 600,
         device = "tiff", compression = "lzw")
  message("Saved: figures/Figure_9B_KEGG_MYC_Correlated_Hub_Genes.tiff")
} else {
  warning("KEGG enrichment returned no significant pathways – Figure 9B skipped.")
}


################################################################################
## 17. FIGURE 10 – GO ENRICHMENT: ME4 Module Hub Genes
################################################################################

message("[17] Running GO enrichment for ME4 module hub genes …")

ego_ME4 <- enrichGO(
  gene          = hub_ME4$GeneSymbol,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_ME4_simp <- clusterProfiler::simplify(ego_ME4, cutoff = 0.7,
                                           by = "p.adjust", select_fun = min)

save_go_dotplot(
  ego_ME4_simp,
  filename = "figures/Figure_10_GO_BP_ME4_Module_Hub_Genes.tiff",
  title    = "GO Biological Processes – MYC-Associated ME4 Module"
)

write.csv(as.data.frame(ego_ME4_simp),
          "tables/06_GO_ME4_module_simplified.csv", row.names = FALSE)


################################################################################
## 18. FIGURES 11A, 11B & 12 – GSEA: MYC Targets + Hallmark Pathways
#
#  pvalueCutoff = 1 retains ALL gene sets so gseaplot2() can always locate
#  V1 / V2 for Figures 11A & 11B.  Figure 12 applies p.adjust < 0.05 filter.
################################################################################

message("[18] Running GSEA …")

gene_list <- res_myc_df %>%
  filter(!is.na(log2FoldChange)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(log2FoldChange, name = gene_symbol)

msig_H <- msigdbr(species = "Homo sapiens", category = "H")

## GSEA for MYC Targets V1 & V2 only
myc_sets <- msig_H %>%
  filter(gs_name %in% c("HALLMARK_MYC_TARGETS_V1",
                         "HALLMARK_MYC_TARGETS_V2")) %>%
  select(gs_name, gene_symbol)

gsea_myc <- GSEA(geneList     = gene_list,
                  TERM2GENE    = myc_sets,
                  pvalueCutoff = 1,
                  verbose      = FALSE)

## GSEA for extended Hallmark gene sets
hallmark_ext <- msig_H %>%
  filter(gs_name %in% c("HALLMARK_MYC_TARGETS_V1",
                         "HALLMARK_MYC_TARGETS_V2",
                         "HALLMARK_E2F_TARGETS",
                         "HALLMARK_G2M_CHECKPOINT",
                         "HALLMARK_DNA_REPAIR")) %>%
  select(gs_name, gene_symbol)

gsea_hallmark <- GSEA(geneList     = gene_list,
                       TERM2GENE    = hallmark_ext,
                       pvalueCutoff = 1,
                       verbose      = FALSE)

## ── Figure 11A: MYC Targets V1 enrichment plot ───────────────────────────────
if ("HALLMARK_MYC_TARGETS_V1" %in% gsea_myc@result$ID) {
  p_gsea_v1 <- gseaplot2(
    gsea_myc,
    geneSetID = "HALLMARK_MYC_TARGETS_V1",
    title     = "GSEA – MYC Targets V1 (CRC)",
    base_size = 13,
    color     = "#D73027"
  )
  tiff("figures/Figure_11A_GSEA_MYC_Targets_V1.tiff",
       width = 7, height = 5, units = "in", res = 600, compression = "lzw")
  print(p_gsea_v1)
  dev.off()
  message("Saved: figures/Figure_11A_GSEA_MYC_Targets_V1.tiff")
} else {
  warning("HALLMARK_MYC_TARGETS_V1 not in GSEA results – Figure 11A skipped.")
}

## ── Figure 11B: MYC Targets V2 enrichment plot ───────────────────────────────
if ("HALLMARK_MYC_TARGETS_V2" %in% gsea_myc@result$ID) {
  p_gsea_v2 <- gseaplot2(
    gsea_myc,
    geneSetID = "HALLMARK_MYC_TARGETS_V2",
    title     = "GSEA – MYC Targets V2 (CRC)",
    base_size = 13,
    color     = "#D73027"
  )
  tiff("figures/Figure_11B_GSEA_MYC_Targets_V2.tiff",
       width = 7, height = 5, units = "in", res = 600, compression = "lzw")
  print(p_gsea_v2)
  dev.off()
  message("Saved: figures/Figure_11B_GSEA_MYC_Targets_V2.tiff")
} else {
  warning("HALLMARK_MYC_TARGETS_V2 not in GSEA results – Figure 11B skipped.")
}

## ── Figure 12: Summary dotplot – significant pathways only ───────────────────
gsea_hallmark_sig        <- gsea_hallmark
gsea_hallmark_sig@result <- gsea_hallmark@result %>% filter(p.adjust < 0.05)

if (nrow(gsea_hallmark_sig@result) > 0) {
  p_gsea_dot <- dotplot(gsea_hallmark_sig, showCategory = 5) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(10, 30, 10, 10)) +
    ggtitle("GSEA – MYC-Associated Hallmark Pathways (CRC)")
  ggsave("figures/Figure_12_GSEA_Hallmark_Pathways_Dotplot.tiff",
         p_gsea_dot, dpi = 600, width = 7, height = 5,
         device = "tiff", compression = "lzw")
  message("Saved: figures/Figure_12_GSEA_Hallmark_Pathways_Dotplot.tiff")
} else {
  warning("No significant hallmark pathways at p.adjust < 0.05 – Figure 12 skipped.")
}

write.csv(as.data.frame(gsea_hallmark),
          "tables/07_GSEA_Hallmark_results.csv", row.names = FALSE)


################################################################################
## 19. FIGURE 13 – SURVIVAL: ME4 Score Kaplan–Meier + Cox Regression
################################################################################

message("[19] Running survival analysis (ME4 score) …")

clinical_surv <- coldata[tumor_mask, ] %>%
  mutate(
    patient_id = rownames(.),
    ME4_score  = ME4_eigengene[patient_id],
    OS_time    = ifelse(!is.na(days_to_death),
                        days_to_death, days_to_last_follow_up),
    OS_event   = ifelse(vital_status == "Dead", 1L, 0L)
  ) %>%
  filter(!is.na(OS_time), OS_time > 0)

message(sprintf("[19] Survival cohort: %d samples  events: %d",
                nrow(clinical_surv),
                sum(clinical_surv$OS_event, na.rm = TRUE)))

clinical_tumor_surv <- clinical_surv %>%
  filter(!is.na(ME4_score))

stopifnot(
  "Insufficient tumor patients for ME4 survival analysis" =
    nrow(clinical_tumor_surv) > 10
)

clinical_tumor_surv <- clinical_tumor_surv %>%
  mutate(ME4_group = factor(
    ifelse(ME4_score >= median(ME4_score, na.rm = TRUE),
           "ME4_high", "ME4_low"),
    levels = c("ME4_low", "ME4_high")
  ))

fit_ME4 <- survfit(Surv(OS_time, OS_event) ~ ME4_group,
                    data = clinical_tumor_surv)

save_tiff_km(
  fit          = fit_ME4,
  data         = clinical_tumor_surv,
  filename     = "figures/Figure_13_KaplanMeier_ME4_Overall_Survival.tiff",
  title        = "Overall Survival – ME4 High vs Low (CRC)",
  legend_title = "ME4 Score",
  legend_labs  = c("ME4 Low", "ME4 High"),
  palette      = c("#4575B4", "#D73027")
)

## Continuous Cox model
cox_ME4_cont <- coxph(Surv(OS_time, OS_event) ~ ME4_score,
                       data = clinical_tumor_surv)
message("[19] Cox PH (continuous ME4 score):")
print(summary(cox_ME4_cont))

## Stage-adjusted Cox model
clinical_tumor_surv <- clinical_tumor_surv %>%
  mutate(stage_simple = factor(substr(ajcc_pathologic_stage, 1, 7)))

df_stage <- clinical_tumor_surv %>%
  filter(!is.na(stage_simple), stage_simple != "")

if (nrow(df_stage) > 20 &&
    nlevels(droplevels(df_stage$stage_simple)) > 1) {
  cox_ME4_stage <- coxph(
    Surv(OS_time, OS_event) ~ ME4_score + stage_simple,
    data = df_stage
  )
  message("[19] Cox PH (ME4 score + AJCC stage):")
  print(summary(cox_ME4_stage))
} else {
  warning("Insufficient or single-level stage data – skipping adjusted Cox.")
}

write.csv(
  clinical_tumor_surv %>%
    select(patient_id, MYC_group, ME4_score, ME4_group,
           OS_time, OS_event, stage_simple),
  "tables/08_Survival_ME4_clinical.csv",
  row.names = FALSE
)


################################################################################
## 20. FIGURE 14 – ML DATASET: Leakage-Free Train/Test Split
################################################################################

message("[20] Building ML dataset with leakage-free 70/30 split …")

## ── Helper: build ML data frame ──────────────────────────────────────────────
build_ml_df <- function(expr_mat, sample_idx) {
  df           <- as.data.frame(t(expr_mat))
  df$MYC_group <- as.character(coldata$MYC_group[sample_idx])
  df$MYC_label <- ifelse(df$MYC_group == "MYC_high", 1L, 0L)
  df$OS_time   <- ifelse(!is.na(coldata$days_to_death[sample_idx]),
                          coldata$days_to_death[sample_idx],
                          coldata$days_to_last_follow_up[sample_idx])
  df$OS_event  <- ifelse(coldata$vital_status[sample_idx] == "Dead", 1L, 0L)
  df
}

## ── A. Full cohort export ─────────────────────────────────────────────────────
tumor_idx_int  <- which(coldata$condition == "Tumor")
expr_full_ME4  <- ensembl_to_symbols(norm_expr[ME4_genes, tumor_idx_int])

full_df <- build_ml_df(expr_full_ME4, tumor_idx_int)
full_df <- full_df[!is.na(full_df$OS_time) & full_df$OS_time > 0, ]

message(sprintf("[20] Full ME4 ML dataset: %d samples × %d genes",
                nrow(full_df), ncol(full_df) - 4))

write.csv(full_df, "ML_dataset_ME4_GeneSymbols_Survival.csv",
          row.names = TRUE)
message("Exported: ML_dataset_ME4_GeneSymbols_Survival.csv")

## ── B. Leakage-free stratified 70/30 split ───────────────────────────────────
set.seed(123)
train_idx <- unlist(lapply(
  split(tumor_idx_int, coldata$MYC_group[tumor_idx_int]),
  function(idx) sample(idx, floor(0.7 * length(idx)))
))
test_idx  <- setdiff(tumor_idx_int, train_idx)

message(sprintf("[20] ML split – train: %d  test: %d",
                length(train_idx), length(test_idx)))

## Run WGCNA on training set only (leakage-free)
expr_train_raw <- norm_expr[, train_idx]
mad_train      <- apply(expr_train_raw, 1, mad)
top5k_train    <- names(sort(mad_train, decreasing = TRUE))[1:5000]
datExpr_train  <- t(expr_train_raw[top5k_train, ])

cor <- WGCNA::cor       # namespace swap IN
net_train <- blockwiseModules(
  datExpr_train,
  power             = softPower,
  nThreads          = MAX_THREADS,   # 6-thread cap
  TOMType           = "unsigned",
  minModuleSize     = 50,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 0
)
cor <- stats::cor       # namespace swap OUT

moduleColors_train <- labels2colors(net_train$colors)
MEs_train          <- moduleEigengenes(datExpr_train,
                                        moduleColors_train)$eigengenes
MYC_num_train      <- as.numeric(coldata$MYC_group[train_idx] == "MYC_high")

modCor_train    <- WGCNA::cor(MEs_train, MYC_num_train, use = "p")
mod_train       <- rownames(modCor_train)[which.max(abs(modCor_train[, 1]))]
color_train     <- gsub("^ME", "", mod_train)
ME4_genes_train <- colnames(datExpr_train)[moduleColors_train == color_train]

message(sprintf("[20] Train ME4 module (%s): %d genes",
                mod_train, length(ME4_genes_train)))

expr_train_ME4 <- ensembl_to_symbols(norm_expr[ME4_genes_train, train_idx])
expr_test_ME4  <- ensembl_to_symbols(norm_expr[ME4_genes_train, test_idx])
common_syms    <- intersect(rownames(expr_train_ME4), rownames(expr_test_ME4))
expr_train_ME4 <- expr_train_ME4[common_syms, ]
expr_test_ME4  <- expr_test_ME4[common_syms, ]

message(sprintf("[20] Shared gene symbols after annotation: %d",
                length(common_syms)))

train_df <- build_ml_df(expr_train_ME4, train_idx)
test_df  <- build_ml_df(expr_test_ME4,  test_idx)

message(sprintf("[20] Train: %d × %d  |  Test: %d × %d",
                nrow(train_df), ncol(train_df) - 4,
                nrow(test_df),  ncol(test_df)  - 4))

write.csv(train_df, "tables/09_ML_train_ME4.csv", row.names = TRUE)
write.csv(test_df,  "tables/10_ML_test_ME4.csv",  row.names = TRUE)

## ── Figure 14: Class balance barplot ─────────────────────────────────────────
bal_df <- bind_rows(
  data.frame(Split = "Train", MYC_group = train_df$MYC_group),
  data.frame(Split = "Test",  MYC_group = test_df$MYC_group)
)
p_bal <- ggplot(bal_df, aes(Split, fill = MYC_group)) +
  geom_bar(position = "fill", width = 0.5) +
  scale_fill_manual(values = c(MYC_low = "#4575B4", MYC_high = "#D73027")) +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "ML Split – Class Balance (ME4 Genes)",
       y = "Proportion", x = NULL, fill = "MYC Group") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/Figure_14_ML_TrainTest_Class_Balance.tiff",
       p_bal, dpi = 600, width = 5, height = 4,
       device = "tiff", compression = "lzw")
message("Saved: figures/Figure_14_ML_TrainTest_Class_Balance.tiff")


################################################################################
## 21. FIGURES 15 & 16 – IMMUNE INFILTRATION: ssGSEA
################################################################################

message("[21] Running ssGSEA immune infiltration analysis …")

immune_sets <- list(
  CD8_T_cells     = c("CD8A", "CD8B"),
  CD4_T_cells     = c("CD4"),
  NK_cells        = c("KLRD1", "NCR1", "NKG7"),
  B_cells         = c("CD19", "MS4A1"),
  Macrophages     = c("CD68", "CSF1R"),
  Dendritic_cells = c("ITGAX", "HLA-DRA"),
  Tregs           = c("FOXP3", "IL2RA"),
  Neutrophils     = c("S100A8", "S100A9")
)

expr_symbols <- ensembl_to_symbols(norm_expr)

ssgsea_param  <- ssgseaParam(
  exprData = as.matrix(expr_symbols),
  geneSets = immune_sets
)
ssgsea_scores <- gsva(ssgsea_param)
immune_scores <- as.data.frame(t(ssgsea_scores))

immune_scores$MYC_group <- coldata$MYC_group[
  match(rownames(immune_scores), rownames(coldata))
]

## ── Figure 15: Immune infiltration heatmap ───────────────────────────────────
immune_mat   <- t(as.matrix(immune_scores[, 1:8]))
sample_order <- order(immune_scores$MYC_group)
immune_mat   <- immune_mat[, sample_order]

annotation_col <- data.frame(
  MYC_group = immune_scores$MYC_group[sample_order],
  row.names  = colnames(immune_mat)
)
annotation_colors <- list(
  MYC_group = c(MYC_low = "#1B9E77", MYC_high = "#E6AB02")
)

tiff("figures/Figure_15_Immune_Infiltration_Heatmap.tiff",
     width = 9, height = 5, units = "in", res = 600, compression = "lzw")
pheatmap(
  scale(immune_mat),
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  main              = "Immune Cell Infiltration – CRC Tumors (MYC Stratified)"
)
dev.off()
message("Saved: figures/Figure_15_Immune_Infiltration_Heatmap.tiff")

## ── Figure 16: Immune score boxplots by MYC group ───────────────────────────
immune_long <- immune_scores %>%
  pivot_longer(cols = -MYC_group,
               names_to  = "Immune_cell",
               values_to = "Score")

p_box <- ggplot(immune_long,
                aes(x = MYC_group, y = Score, fill = MYC_group)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~Immune_cell, scales = "free") +
  scale_fill_manual(values = c(MYC_low = "#1B9E77", MYC_high = "#E6AB02")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Immune Infiltration Differences: MYC-high vs MYC-low",
    x     = "MYC group",
    y     = "ssGSEA score"
  ) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90")
  )

ggsave("figures/Figure_16_Immune_MYCHigh_vs_Low_Boxplots.tiff",
       p_box, dpi = 600, width = 9, height = 6,
       device = "tiff", compression = "lzw")
message("Saved: figures/Figure_16_Immune_MYCHigh_vs_Low_Boxplots.tiff")


################################################################################
## 22. FIGURE 17 – MODULE vs IMMUNE INFILTRATION
################################################################################

message("[22] Correlating WGCNA turquoise module with immune infiltration …")

module_scores  <- MEs$MEturquoise
names(module_scores) <- rownames(MEs)

immune_samples <- substr(rownames(immune_scores), 1, 12)
module_samples <- substr(names(module_scores),    1, 12)

module_aligned <- module_scores[match(immune_samples, module_samples)]
immune_scores$Turquoise <- module_aligned

immune_matrix <- immune_scores[, c(
  "CD8_T_cells", "CD4_T_cells", "NK_cells", "B_cells",
  "Macrophages", "Dendritic_cells", "Tregs", "Neutrophils"
)]

cor_results <- apply(immune_matrix, 2, function(x)
  cor(x, immune_scores$Turquoise,
      method = "spearman",
      use    = "pairwise.complete.obs"))

cor_mat            <- as.matrix(cor_results)
colnames(cor_mat)  <- "Turquoise_module"

tiff("figures/Figure_17_Turquoise_Immune_Correlation.tiff",
     width = 7, height = 7, units = "in", res = 600, compression = "lzw")
pheatmap(
  cor_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color        = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  main         = "MYC Module vs Immune Infiltration",
  fontsize     = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col    = "0",
  border_color = NA
)
dev.off()
message("Saved: figures/Figure_17_Turquoise_Immune_Correlation.tiff")


################################################################################
## 23. FIGURE 18 – IMMUNE CHECKPOINT EXPRESSION: MYC-High vs MYC-Low
################################################################################

message("[23] Plotting immune checkpoint expression …")

checkpoint_genes <- c(
  "PDCD1",   # PD-1
  "CD274",   # PD-L1
  "CTLA4",
  "LAG3",
  "TIGIT",
  "HAVCR2"   # TIM-3
)

checkpoint_expr <- expr_symbols[
  rownames(expr_symbols) %in% checkpoint_genes, , drop = FALSE
]

checkpoint_df <- as.data.frame(t(checkpoint_expr))
checkpoint_df$MYC_group <- coldata$MYC_group[
  match(rownames(checkpoint_df), rownames(coldata))
]

checkpoint_long <- checkpoint_df %>%
  pivot_longer(cols = -MYC_group,
               names_to  = "Checkpoint",
               values_to = "Expression")

p_checkpoint <- ggplot(checkpoint_long,
                        aes(x = MYC_group, y = Expression, fill = MYC_group)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~Checkpoint, scales = "free") +
  scale_fill_manual(values = c(MYC_low = "#1B9E77", MYC_high = "#E6AB02")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Immune Checkpoint Expression: MYC-high vs MYC-low CRC",
    x     = "MYC group",
    y     = "Expression"
  ) +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90")
  )

ggsave("figures/Figure_18_Immune_Checkpoints_MYC.tiff",
       p_checkpoint, dpi = 600, width = 9, height = 6,
       device = "tiff", compression = "lzw")
message("Saved: figures/Figure_18_Immune_Checkpoints_MYC.tiff")


################################################################################
## 24. FIGURE 19 – MYC–IMMUNE INTERACTION NETWORK
################################################################################

message("[24] Building MYC–Immune interaction network …")

immune_cells <- c("CD8_T_cells", "CD4_T_cells", "NK_cells", "B_cells",
                   "Macrophages", "Dendritic_cells", "Tregs", "Neutrophils")

## Compute correlations for all immune + checkpoint features
all_features <- c(immune_cells, checkpoint_genes)

cor_df <- data.frame(
  Feature     = all_features,
  Correlation = sapply(all_features, function(feat) {
    if (feat %in% colnames(immune_scores)) {
      cor(immune_scores[[feat]], immune_scores$Turquoise,
          method = "spearman", use = "pairwise.complete.obs")
    } else if (feat %in% rownames(checkpoint_expr)) {
      cor(as.numeric(checkpoint_expr[feat, ]),
          immune_scores$Turquoise[match(colnames(checkpoint_expr),
                                         rownames(immune_scores))],
          method = "spearman", use = "pairwise.complete.obs")
    } else {
      NA_real_
    }
  })
) %>%
  filter(!is.na(Correlation))

cor_df$type <- ifelse(cor_df$Feature %in% immune_cells,
                       "Immune_cell", "Checkpoint")

edges <- cor_df %>%
  filter(abs(Correlation) > 0.2) %>%
  transmute(
    from        = "MYC_module",
    to          = Feature,
    Correlation = Correlation,
    sign        = ifelse(Correlation > 0, "Positive", "Negative")
  )

nodes <- data.frame(
  name = unique(c("MYC_module", edges$to)),
  stringsAsFactors = FALSE
) %>%
  mutate(type = case_when(
    name == "MYC_module"        ~ "MYC_module",
    name %in% immune_cells      ~ "Immune_cell",
    TRUE                        ~ "Checkpoint"
  ))

g <- graph_from_data_frame(
  edges[, c("from", "to", "Correlation", "sign")],
  vertices = nodes
)

p_net <- ggraph(g, layout = "fr") +
  geom_edge_link(
    aes(color = sign, width = abs(Correlation)),
    alpha = 0.8
  ) +
  geom_node_point(aes(color = type), size = 7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_color_manual(
    values = c(
      MYC_module  = "#D73027",
      Immune_cell = "#4575B4",
      Checkpoint  = "#66A61E",
      Positive    = "#D73027",
      Negative    = "#4575B4"
    )
  ) +
  scale_edge_color_manual(
    values = c(Positive = "#D73027", Negative = "#4575B4")
  ) +
  scale_edge_width(range = c(0.5, 2.5)) +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white"),
    legend.position  = "right"
  ) +
  ggtitle("MYC Module–Immune Interaction Network in CRC")

ggsave("figures/Figure_19_MYC_Immune_Network.tiff",
       p_net, width = 7, height = 6, dpi = 600,
       device = "tiff", compression = "lzw")
message("Saved: figures/Figure_19_MYC_Immune_Network.tiff")


################################################################################
## 25. TABLES 11 & 12 – ceRNA NETWORK  (miRTarBase + ENCORI)
#
#  Requires local copies of:
#    miRTarBase_MTI.csv       https://mirtarbase.cuhk.edu.cn
#    ENCORI_miRNA_lncRNA.tsv  https://rnasysu.com/encori
################################################################################

message("[25] Building ceRNA network …")

MIRTARBASE_FILE <- "miRTarBase_MTI.csv"
ENCORI_FILE     <- "ENCORI_miRNA_lncRNA.tsv"

if (!file.exists(MIRTARBASE_FILE)) {
  warning(
    "miRTarBase file not found – ceRNA network skipped.\n",
    "Download from: https://mirtarbase.cuhk.edu.cn"
  )
} else if (!file.exists(ENCORI_FILE)) {
  warning(
    "ENCORI file not found – ceRNA network skipped.\n",
    "Download from: https://rnasysu.com/encori"
  )
} else {

  myc_mRNAs <- res_myc_df %>%
    filter(padj < 0.05, log2FoldChange > 1, !is.na(gene_symbol)) %>%
    distinct(gene_symbol) %>%
    slice_head(n = 40) %>%
    pull(gene_symbol)

  mirtar   <- read.csv(MIRTARBASE_FILE, stringsAsFactors = FALSE)
  mir_mrna <- mirtar %>%
    dplyr::filter(grepl("^hsa", miRNA),
                  Target.Gene   %in% myc_mRNAs,
                  Support.Type  == "Functional MTI") %>%
    dplyr::select(miRNA, mRNA = Target.Gene) %>%
    dplyr::distinct()

  message(sprintf("[25] miRNA–mRNA interactions: %d", nrow(mir_mrna)))

  mir_lnc_raw   <- fread(ENCORI_FILE, data.table = FALSE)
  mir_lnc_clean <- mir_lnc_raw %>%
    dplyr::filter(geneType   == "lncRNA",
                  miRNAname  %in% mir_mrna$miRNA) %>%
    dplyr::select(lncRNA = geneName, miRNA = miRNAname) %>%
    dplyr::distinct()

  message(sprintf("[25] lncRNA–miRNA interactions: %d", nrow(mir_lnc_clean)))

  ceRNA_triplets <- mir_lnc_clean %>%
    inner_join(mir_mrna, by = "miRNA") %>%
    distinct()

  message(sprintf("[25] ceRNA triplets: %d", nrow(ceRNA_triplets)))

  edges_cerna <- bind_rows(
    ceRNA_triplets %>%
      transmute(source = lncRNA, target = miRNA,
                interaction = "lncRNA-miRNA"),
    ceRNA_triplets %>%
      transmute(source = miRNA,  target = mRNA,
                interaction = "miRNA-mRNA")
  )

  nodes_cerna <- data.frame(
    id = unique(c(edges_cerna$source, edges_cerna$target)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(type = case_when(
      grepl("^hsa",                              id) ~ "miRNA",
      grepl("^LINC|^MALAT|^NEAT|^HOTAIR|^OIP5", id) ~ "lncRNA",
      TRUE                                           ~ "mRNA"
    ))

  write.csv(edges_cerna, "tables/11_ceRNA_network_edges.csv", row.names = FALSE)
  write.csv(nodes_cerna, "tables/12_ceRNA_network_nodes.csv", row.names = FALSE)

  message(sprintf("[25] ceRNA network: %d edges, %d nodes",
                  nrow(edges_cerna), nrow(nodes_cerna)))
}


################################################################################
message("")
message("══════════════════════════════════════════════════════════════════════")
message("  Pipeline complete.")
message("  Figures : figures/  (Figures 1–19)")
message("  Tables  : tables/   (Tables 01–12)")
message("══════════════════════════════════════════════════════════════════════")
################################################################################
