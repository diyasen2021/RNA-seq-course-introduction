# GEO RNA-seq Pipeline
## Part 4 — DESeq2 Differential Expression from a GEO Count Matrix

> **Note on continuity:** Parts 1–3 walked through downloading FASTQs, running FastQC/Trimmomatic, and quantifying with Salmon yourself. Part 4 makes a deliberate detour: instead of continuing from your own `quant.sf` files, we go back to GEO and pull NCBI's pre-built gene-level count matrix (`GSEnnnnn_raw_counts_GRCh38.p13_NCBI.tsv.gz`) as the input. This is a completely valid and common way to start a DESeq2 analysis — NCBI re-processes the raw SRA reads for every GEO RNA-seq series through a standardized STAR + RSEM pipeline, so what you get back is already integer, gene-level (Entrez ID) raw counts. No tximeta, no transcript-to-gene collapsing needed — we go straight from matrix to DESeq2.

## Contents

1. [Why start from the GEO matrix instead of your own Salmon output](#1-why-start-from-the-geo-matrix-instead-of-your-own-salmon-output)
2. [Understanding the NCBI-generated count matrix](#2-understanding-the-ncbi-generated-count-matrix)
3. [Setting up your R environment](#3-setting-up-your-r-environment)
4. [Loading the count matrix](#4-loading-the-count-matrix)
5. [Building sample metadata — the paired design](#5-building-sample-metadata--the-paired-design)
6. [Why this is a paired (repeated-measures) design, not an unpaired one](#6-why-this-is-a-paired-repeated-measures-design-not-an-unpaired-one)
7. [Constructing the DESeqDataSet](#7-constructing-the-deseqdataset)
8. [Pre-filtering low-count genes](#8-pre-filtering-low-count-genes)
9. [Running DESeq2](#9-running-deseq2)
10. [Quality control: PCA and sample distances](#10-quality-control-pca-and-sample-distances)
11. [Extracting pairwise contrasts](#11-extracting-pairwise-contrasts)
12. [Volcano plots for each contrast](#12-volcano-plots-for-each-contrast)
13. [Saving ranked DEG tables](#13-saving-ranked-deg-tables)
14. [Part 4 checkpoint](#14-part-4-checkpoint)

---

## 1. Why start from the GEO matrix instead of your own Salmon output

There are two legitimate entry points into a DESeq2 analysis, and it's worth being explicit about which one we're using and why.

**Entry point A — your own Salmon `quant.sf` files (Parts 1–3):** you control every step from raw FASTQ to quantification. This is the right approach when you need to justify and document your own QC/trimming/alignment decisions, or when the GEO submitter didn't provide a usable raw count matrix.

**Entry point B — NCBI's re-processed count matrix (this lesson):** NCBI downloads the original SRA reads for every RNA-seq GEO series and re-runs them through one standardized pipeline (STAR alignment to GRCh38.p13 + RSEM quantification), publishing the result as `GSEnnnnn_raw_counts_GRCh38.p13_NCBI.tsv.gz` on the GEO accession page. This is useful when you want a fast, reproducible starting point, when you want to sanity-check your own Salmon counts against an independent pipeline, or when you're focusing your effort on the downstream statistics and biology rather than re-deriving counts that already exist.

Both are "correct" — they are just different, equally valid pipelines (pseudoalignment+EM vs. alignment+RSEM) converging on the same kind of answer: an integer count matrix suitable for DESeq2. For this lesson we use Entry point B. If you also have Salmon counts from Part 3, comparing the two is a genuinely good portfolio exercise — but that is optional and not required to proceed.

---

## 2. Understanding the NCBI-generated count matrix

On the GEO accession page, scroll to the **"Supplementary file"** section. You're looking for a filename matching this pattern:

```
GSEnnnnn_raw_counts_GRCh38.p13_NCBI.tsv.gz
```

Download and inspect it before loading into R:

```bash
gunzip -k GSEnnnnn_raw_counts_GRCh38.p13_NCBI.tsv.gz
head -5 GSEnnnnn_raw_counts_GRCh38.p13_NCBI.tsv
```

```
GeneID    GSMxxxxxx1  GSMxxxxxx2  GSMxxxxxx3  GSMxxxxxx4  ...
100287102 0           0           0           0
100287596 12          8           15          3
100288569 0           1           0           0
653635    142         98          201         87
```

A few things to note about this file before you touch it in R:

- **Rows are genes, identified by Entrez Gene ID** — not Ensembl IDs, and not gene symbols. You will need to convert these later (e.g. with `org.Hs.eg.db` or `AnnotationDbi`) before the results are human-readable.
- **Columns are samples, identified by GSM accession** — the individual sample ID assigned by GEO, *not* your own sample naming scheme. You will map these to your own sample metadata in Section 5.
- **Values are raw integer counts** — this is exactly what DESeq2 wants. No TPM, no normalization, no decimals. Confirm this yourself the first time you use a new GEO matrix: `class(counts_matrix[1,1])` should show integer-like values, and `any(counts_matrix != round(counts_matrix))` should be `FALSE`.
- **The number of rows is large (~40,000+)** — this includes many genes with zero or near-zero counts across all samples. Section 8 covers filtering these out before running DESeq2.

---

## 3. Setting up your R environment

```r
# Bioconductor packages (install once, if not already present)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", "apeglm"))

# CRAN packages for plotting
install.packages(c("ggplot2", "ggrepel", "pheatmap"))

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
```

`apeglm` is included here because it gives better (less noisy) log fold-change shrinkage than the default DESeq2 estimator, which matters when you have a relatively small number of samples per group — exactly the situation in a typical GEO cohort.

---

## 4. Loading the count matrix

```r
# Read the NCBI-generated count matrix
# header = TRUE because row 1 contains GSM sample IDs
# row.names = 1 because column 1 is the Entrez GeneID
counts_raw <- read.delim(
  "data/raw_counts/GSEnnnnn_raw_counts_GRCh38.p13_NCBI.tsv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE   # preserves GSM IDs exactly as written (avoids R mangling "-" etc.)
)

# Convert to a plain matrix of integers — DESeq2 requires this, not a data.frame
counts_matrix <- as.matrix(counts_raw)

# Sanity checks — always run these before proceeding
dim(counts_matrix)                                  # genes x samples
class(counts_matrix[1, 1])                           # should be numeric/integer
all(counts_matrix == round(counts_matrix))            # should be TRUE — confirms raw integer counts
colnames(counts_matrix)                               # confirm these are GSM accessions
```

If `all(counts_matrix == round(counts_matrix))` returns `FALSE`, stop here — you have been given a normalized or otherwise transformed matrix, not raw counts, and DESeq2 will produce silently incorrect results if you proceed.

---

## 5. Building sample metadata — the paired design

The count matrix only tells you *what* was measured. It says nothing about which sample belongs to which patient or which condition — that information lives in the GEO **series matrix** / sample characteristics, and you have to build it yourself as a metadata table.

This is written generically below — using a small number of placeholder variables you set once at the top — so you can adapt it to any cohort size without rewriting the logic.

```r
# --- Set these once, based on your specific GEO series ---
n_patients        <- length(unique(condition_levels))   # placeholder, replaced below
condition_levels  <- c("normal", "primary", "metastasis")  # edit to match your study design
samples_per_patient <- length(condition_levels)

# --- Build metadata programmatically rather than typing it by hand ---
# gsm_ids should be in the SAME order as colnames(counts_matrix)
gsm_ids <- colnames(counts_matrix)

# patient_id and condition vectors must be supplied from the GEO sample
# characteristics (series matrix) — example construction shown below assumes
# samples are grouped in blocks of `samples_per_patient`, in matrix column order.
# Replace this block with an explicit mapping if your matrix columns are not
# already grouped by patient.

n_total_samples <- ncol(counts_matrix)
n_patients      <- n_total_samples / samples_per_patient

coldata <- data.frame(
  gsm_id      = gsm_ids,
  patient_id  = factor(rep(paste0("patient_", seq_len(n_patients)),
                            each = samples_per_patient)),
  condition   = factor(rep(condition_levels, times = n_patients),
                        levels = condition_levels),
  row.names   = gsm_ids
)

head(coldata)
table(coldata$condition)        # should be balanced: same n per condition
table(coldata$patient_id)       # should show samples_per_patient for every patient
```

**Important:** the placeholder construction above assumes your matrix columns are already ordered patient-by-patient, condition-by-condition. In practice, GEO column order is not guaranteed to follow this pattern — always cross-check `coldata` against the actual GEO sample characteristics table (visible on the GEO accession page under each GSM) before proceeding. A mismatched condition/patient label is the single most common cause of meaningless DESeq2 results, and it produces no error — only wrong answers.

```r
# Critical check: do sample names match between matrix and metadata, IN ORDER?
all(colnames(counts_matrix) == rownames(coldata))   # MUST be TRUE before continuing
```

If this returns `FALSE`, do not proceed. Either reorder `coldata` to match `counts_matrix`, or reorder the matrix columns — `counts_matrix <- counts_matrix[, rownames(coldata)]` is the usual fix.

---

## 6. Why this is a paired (repeated-measures) design, not an unpaired one

This point is worth pausing on because it changes the statistics, not just the code.

A naive approach would treat every sample as independent and just compare "all normal samples" vs. "all primary samples" vs. "all metastasis samples" as three unrelated groups. That is wrong here, because each patient contributes one sample to *every* condition. Samples from the same patient share that patient's genetic background, age, comorbidities, RNA extraction batch, and baseline transcriptional state — they are correlated with each other in ways that have nothing to do with the condition being tested.

If you ignore this correlation, patient-to-patient variability gets absorbed into what looks like condition-driven variability, inflating your false discovery rate. The fix is to include `patient_id` as a blocking term in the design formula, which tells DESeq2: "first remove each patient's overall baseline expression level, then test for the effect of condition within each patient."

This is directly analogous to a paired t-test vs. an unpaired t-test — same logic, extended to a count-based GLM framework.

---

## 7. Constructing the DESeqDataSet

```r
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = coldata,
  design    = ~ patient_id + condition
)

dds
```

The design formula `~ patient_id + condition` does two things: `patient_id` absorbs (blocks out) each patient's baseline expression differences, and `condition` — listed last, which matters in DESeq2 syntax — is the term whose effect you actually want to test. DESeq2 always tests the last term in the formula by default.

---

## 8. Pre-filtering low-count genes

```r
# Remove genes with essentially no signal across all samples.
# This is not strictly required (DESeq2 handles low counts internally),
# but it speeds up the analysis and reduces multiple-testing burden.
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

dim(dds)   # confirm the gene count has dropped from ~40,000+ to a smaller, real set
```

The threshold here (`>=10` counts in at least 3 samples) is a reasonable generic default — not a strict rule. Adjust the "3" to roughly match the smallest group size in your `condition_levels` once you know it.

---

## 9. Running DESeq2

```r
dds <- DESeq(dds)

# Inspect the fitted model
resultsNames(dds)
```

`DESeq()` runs three steps internally in one call: estimation of size factors (sample-to-sample normalization, the median-of-ratios method — this is why you must give it raw counts, not TPM), estimation of per-gene dispersion (how much more variable each gene's counts are than a simple Poisson model would predict), and fitting of a negative binomial GLM per gene followed by Wald testing for the `condition` term.

`resultsNames(dds)` will show you the exact coefficient names DESeq2 created from your `condition` factor — you need these exact strings for Section 11.

---

## 10. Quality control: PCA and sample distances

Always look at this before trusting any DE results below.

```r
# Variance-stabilizing transformation — for visualization only, not for DE testing
vsd <- vst(dds, blind = TRUE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = c("condition", "patient_id"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = patient_id)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA of all samples, coloured by condition")
```

What to look for: samples from the same condition should cluster together more than samples from the same patient do — if instead samples cluster strongly by patient regardless of condition, that's exactly why you blocked on `patient_id` in the design formula. Some patient-driven clustering is normal and expected; condition should still be a visible, separate axis of variation (usually PC1 or PC2).

```r
# Sample-to-sample distance heatmap — a second, complementary view of the same question
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

pheatmap(sample_dist_matrix,
         annotation_col = as.data.frame(colData(dds)[, c("condition", "patient_id")]),
         main = "Sample-to-sample distances")
```

---

## 11. Extracting pairwise contrasts

With three (or more) condition levels, you need explicit pairwise contrasts rather than relying on the default comparison alone. Replace the level names below with whatever you set in `condition_levels` in Section 5.

```r
# Contrast 1: second level vs. first level (e.g., primary vs. normal)
res_c1 <- results(dds, contrast = c("condition", condition_levels[2], condition_levels[1]))

# Contrast 2: third level vs. first level (e.g., metastasis vs. normal)
res_c2 <- results(dds, contrast = c("condition", condition_levels[3], condition_levels[1]))

# Contrast 3: third level vs. second level (e.g., metastasis vs. primary)
res_c3 <- results(dds, contrast = c("condition", condition_levels[3], condition_levels[2]))

# Apply shrinkage to log fold-changes for more stable ranking/visualization
# (coef names come from resultsNames(dds) in Section 9 — adjust indices if needed)
res_c1_shrunk <- lfcShrink(dds, contrast = c("condition", condition_levels[2], condition_levels[1]), res = res_c1, type = "ashr")
```

A quick summary of each contrast tells you how many genes passed significance before you go further:

```r
summary(res_c1)
summary(res_c2)
summary(res_c3)
```

Note: `lfcShrink()` with `type = "ashr"` requires the `ashr` package (`install.packages("ashr")`); if you'd rather stick to `apeglm`, that method requires the contrast to be expressed as a model coefficient rather than an arbitrary contrast vector — check `resultsNames(dds)` and use `coef = "..."` instead of `contrast = ...` in that case.

---

## 12. Volcano plots for each contrast

```r
plot_volcano <- function(res_df, title, fc_threshold = 1, padj_threshold = 0.05) {
  res_df <- as.data.frame(res_df)
  res_df$gene <- rownames(res_df)
  res_df$significant <- with(res_df,
    !is.na(padj) & padj < padj_threshold & abs(log2FoldChange) > fc_threshold
  )

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
    theme_bw() +
    ggtitle(title) +
    xlab("log2 fold change") +
    ylab("-log10 adjusted p-value")
}

plot_volcano(res_c1, paste(condition_levels[2], "vs", condition_levels[1]))
plot_volcano(res_c2, paste(condition_levels[3], "vs", condition_levels[1]))
plot_volcano(res_c3, paste(condition_levels[3], "vs", condition_levels[2]))
```

Treat the thresholds (`fc_threshold = 1`, i.e. 2-fold change, and `padj_threshold = 0.05`) as a starting convention, not a fixed rule — state explicitly in your methods section whatever thresholds you actually used.

---

## 13. Saving ranked DEG tables

```r
# Helper: convert Entrez IDs to gene symbols and write a clean, ranked CSV
save_deg_table <- function(res, filename) {
  res_df <- as.data.frame(res)
  res_df$entrez_id <- rownames(res_df)

  res_df$gene_symbol <- mapIds(
    org.Hs.eg.db,
    keys = res_df$entrez_id,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )

  res_df <- res_df[order(res_df$padj), ]
  res_df <- res_df[, c("entrez_id", "gene_symbol", "baseMean",
                        "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

  write.csv(res_df, filename, row.names = FALSE)
  invisible(res_df)
}

deg_c1 <- save_deg_table(res_c1, "results/deseq2/contrast1_DEGs.csv")
deg_c2 <- save_deg_table(res_c2, "results/deseq2/contrast2_DEGs.csv")
deg_c3 <- save_deg_table(res_c3, "results/deseq2/contrast3_DEGs.csv")

# Quick look at the top genes from the most clinically interesting contrast
head(deg_c3, 15)
```

---

## 14. Part 4 checkpoint

Before moving to Part 5, verify all of the following:

- [ ] Raw count matrix loaded and confirmed to contain integer values (not TPM/FPKM)
- [ ] `colnames(counts_matrix)` matches `rownames(coldata)` exactly, in the same order
- [ ] Sample metadata correctly encodes `patient_id` and `condition` for every sample, cross-checked against the GEO sample characteristics table
- [ ] DESeqDataSet built with design `~ patient_id + condition`
- [ ] Low-count genes filtered before running `DESeq()`
- [ ] PCA plot reviewed — condition forms a visible axis of separation
- [ ] All pairwise contrasts of interest extracted with `results()`
- [ ] Log fold-changes shrunk with `lfcShrink()` before ranking/plotting
- [ ] Volcano plots generated for each contrast
- [ ] Ranked DEG tables saved as CSV with gene symbols (not just Entrez IDs)
- [ ] All thresholds used (padj, log2FC) documented in your methods notes
- [ ] All changes committed to git

```bash
git add scripts/r/04_deseq2_analysis.R results/deseq2/ 
git commit -m "Part 4 complete: DESeq2 differential expression from GEO count matrix"
```

---

## Coming up in Part 5

Part 5 takes the ranked DEG tables produced here into pathway-level interpretation — Hallmark gene set enrichment, and a closer look at the genes specifically distinguishing metastatic progression from primary disease.

| Part | Title | What you'll do |
|------|-------|----------------|
| ✅ Part 1 | Setup & Download | WSL, conda, SRA download, Salmon index |
| ✅ Part 2 | QC & Trimming | FastQC, MultiQC, Trimmomatic |
| ✅ Part 3 | Salmon Quantification | Pseudoalignment, quant.sf files, mapping QC |
| ✅ Part 4 | DESeq2 Differential Expression | GEO count matrix import, paired design, three contrasts, volcano plots |
| **▶ Part 5** | **GSEA & Pathway Analysis** | **Hallmark enrichment, progression gene signatures, figures** |

---

*GEO RNA-seq Pipeline · Part 4 of 5 · Masters in Bioinformatics Portfolio Project*
