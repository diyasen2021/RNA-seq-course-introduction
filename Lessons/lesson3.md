# GSE50760 CRC RNA-seq Pipeline
## Part 3 — Salmon Quantification

## Contents

1. [What is Salmon and why are we using it](#1-what-is-salmon-and-why-are-we-using-it)
2. [How pseudoalignment works](#2-how-pseudoalignment-works)
3. [Understanding the quant.sf output file](#3-understanding-the-quantsf-output-file)
4. [TPM vs raw counts — a critical distinction](#4-tpm-vs-raw-counts--a-critical-distinction)
5. [Library type detection for GSE50760](#5-library-type-detection-for-gse50760)
6. [Running Salmon on all 18 samples](#6-running-salmon-on-all-18-samples)
7. [Checking mapping rates](#7-checking-mapping-rates)
8. [Aggregating Salmon QC with MultiQC](#8-aggregating-salmon-qc-with-multiqc)
9. [Verifying output structure](#9-verifying-output-structure)
10. [Introduction to tximeta — what comes next](#10-introduction-to-tximeta--what-comes-next)
11. [Part 3 checkpoint](#11-part-3-checkpoint)

---

## 1. What is Salmon and why are we using it

At the end of Part 2 you have 18 trimmed FASTQ files — one per sample, covering 6 normal colon, 6 primary CRC, and 6 liver metastasis samples from the AMC cohort. The next question is: **for each gene in the human genome, how many reads came from that gene in each sample?**

There are two broad approaches to answering this:

**Alignment-based:** Map every read to the genome using a splice-aware aligner (STAR or HISAT2) → produce a BAM file → count reads overlapping each gene with featureCounts or HTSeq. This is conceptually intuitive but requires 32 GB RAM for the STAR genome index, produces 5–15 GB BAM files per sample, and takes 30–60 minutes per sample on a laptop.

**Pseudoalignment-based (what we use):** Skip the genome entirely. Map reads directly against the transcriptome (cDNA sequences) using k-mer hashing → produce transcript-level count estimates directly. Salmon is the leading tool for this. It requires ~8 GB RAM, produces tiny output files (<1 MB per sample), and runs in 2–5 minutes per sample.

This is not a shortcut — pseudoalignment is the current best practice for differential expression analysis. It is more accurate than alignment-based counting for multi-mapping reads and produces counts that are statistically better suited to DESeq2. The majority of published RNA-seq DE analyses in 2022 onwards use Salmon or kallisto.

The one thing pseudoalignment cannot do is produce a BAM file for genome browser visualisation. If you need to visually inspect read coverage at specific loci you still need STAR. For this project you do not.

---

## 2. How pseudoalignment works

Understanding the algorithm at the right level of detail is important — you will be asked about this in interviews and thesis defences.

### The k-mer index

When you built the Salmon index in Part 1, Salmon broke every transcript in the Ensembl cDNA reference into overlapping k-mers (short sequences of length k=31 by default) and stored which transcripts each k-mer appears in. This created a hash table: k-mer sequence → set of compatible transcripts.

```
Transcript A: ATCGATCGATCG...
Transcript B: ATCGATCGTTCG...

k-mer "ATCGATCG" → {Transcript A, Transcript B}  (shared)
k-mer "ATCGATTT" → {Transcript A only}             (unique to A)
k-mer "ATCGTTCG" → {Transcript B only}             (unique to B)
```

### Quasi-mapping

For each read in your FASTQ file, Salmon:
1. Extracts k-mers from the read
2. Looks each k-mer up in the index hash table
3. Finds the set of transcripts compatible with the read (the "equivalence class")
4. Does NOT find the exact position within those transcripts — this is what makes it fast

### Expectation-Maximisation (EM)

After all reads are mapped to equivalence classes, Salmon uses an EM algorithm to estimate how many reads came from each individual transcript. This handles multi-mapping reads (reads that could have come from multiple transcripts) probabilistically and correctly — alignment-based methods either discard these reads or assign them naively.

### Why this matters for CRC specifically

The human transcriptome has extensive alternative splicing, and colorectal cancer tumours express cancer-specific splice isoforms. Salmon handles these more accurately than alignment + featureCounts because the EM step properly distributes multi-mapping reads across isoforms rather than discarding them.

---

## 3. Understanding the quant.sf output file

Salmon produces one directory per sample, each containing a `quant.sf` file. Let's look at what that file contains before running anything.

```
quant.sf — tab-delimited, one row per transcript
```

| Column | Name | What it is |
|--------|------|-----------|
| 1 | Name | Ensembl transcript ID (e.g., `ENST00000456328.2`) |
| 2 | Length | Full transcript length in bases |
| 3 | EffectiveLength | Length corrected for fragment size distribution and GC bias |
| 4 | TPM | Transcripts Per Million — length-normalised, depth-normalised |
| 5 | NumReads | **Estimated read count** — what DESeq2 uses |

A real `quant.sf` file looks like this:

```
Name                    Length  EffectiveLength  TPM         NumReads
ENST00000456328.2       1657    1508.21          0.000000    0.000
ENST00000450305.2       632     483.21           0.000000    0.000
ENST00000488147.1       1351    1202.21          12.341280   892.412
ENST00000619216.1       68      0.00             0.000000    0.000
ENST00000473358.1       712     563.21           0.000000    0.000
```

There are ~230,000 rows — one per transcript in the Ensembl reference. When we import into DESeq2 via `tximeta` in Part 4, we will aggregate these transcript-level counts to gene-level counts automatically.

---

## 4. TPM vs raw counts — a critical distinction

This is one of the most common mistakes in RNA-seq analysis. It must be understood before you run a single line of DESeq2.

### TPM (Transcripts Per Million)

TPM normalises for two things:
1. **Gene length** — longer genes have more reads just because they are longer
2. **Sequencing depth** — samples sequenced more deeply have more reads for every gene

TPM is useful for: visualisation, heatmaps, comparing expression levels between genes within a sample, and reporting expression of individual genes.

TPM is **not** suitable for: DESeq2 or edgeR. These tools implement their own normalisation internally (DESeq2 uses median-of-ratios). If you give them TPM values they will normalise already-normalised data and produce statistically wrong results — silently, with no error message.

### NumReads (estimated counts)

These are not integers — they are fractional estimates from the EM algorithm (e.g., `892.412`). They represent the estimated number of reads that came from each transcript. DESeq2 accepts these fractional counts correctly via `tximeta`.

**Rule:** Always use `NumReads` for DESeq2. Always.

### What this means for GSE50760 specifically

The files in `GSE50760_RAW.tar` (the per-sample FPKM files) are **normalised** and **cannot be used for DESeq2**. The NCBI-generated `GSE50760_raw_counts_GRCh38.p13_NCBI.tsv.gz` file contains integer counts and **can** be used directly. The `quant.sf` files you are about to generate using Salmon contain fractional estimated counts and **can** be used via tximeta.

---

## 5. Library type detection for GSE50760

Before running Salmon you need to know whether the GSE50760 libraries are stranded or unstranded. This matters because strandedness affects how Salmon assigns reads that could map to either strand.

From the series matrix file we read in the last session: the library was prepared using **TruSeq RNA Sample Preparation kit v2** with poly-A selection and sequenced **paired-end 2x100bp**. TruSeq v2 is an **unstranded** library prep protocol.

However, rather than hard-coding this, we use `--libType A` which tells Salmon to **auto-detect** the library type from the first 100,000 reads. This is best practice because:
- It provides a sanity check — if Salmon detects something unexpected it warns you
- The auto-detected library type is written to the output and can be verified
- It costs nothing in runtime

Salmon's library type codes:

| Code | Meaning |
|------|---------|
| `A` | Auto-detect (recommended for most cases) |
| `U` | Unstranded |
| `SF` | Stranded forward |
| `SR` | Stranded reverse |
| `ISF` | Paired-end, inward, stranded forward |
| `ISR` | Paired-end, inward, stranded reverse |
| `IU` | Paired-end, inward, unstranded |

For GSE50760 you expect Salmon to report `IU` — paired-end, inward orientation, unstranded. If it reports something different, flag it before proceeding.

---

## 6. Running Salmon on all 18 samples

GSE50760 is **paired-end** — each sample has an `_1.fastq.gz` and `_2.fastq.gz` file. This is different from Part 1 where we discussed single-end reads. The Salmon command changes: instead of `-r` for a single file, we use `-1` and `-2` for the read pairs.

Save this as `scripts/bash/05_salmon_quant.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────
TRIMDIR="data/trimmed"
OUTDIR="results/salmon_quant"
INDEX="reference/salmon_index"
THREADS=6
LOGDIR="logs/salmon"
mkdir -p "$OUTDIR" "$LOGDIR"

# Verify the index exists
if [[ ! -d "$INDEX" ]]; then
  echo "ERROR: Salmon index not found at ${INDEX}"
  echo "Did you complete Part 1? Run: salmon index -t reference/transcriptome/*.fa.gz -i reference/salmon_index -p 6"
  exit 1
fi

echo "[$(date '+%H:%M:%S')] Starting Salmon quantification for all samples"
echo "[$(date '+%H:%M:%S')] Index: ${INDEX}"
echo "[$(date '+%H:%M:%S')] Threads: ${THREADS}"

# ── Sample loop ─────────────────────────────────────────────────────────
# GSE50760 naming: SRR_XXXXXXX_1_trimmed.fastq.gz and SRR_XXXXXXX_2_trimmed.fastq.gz
for R1 in "${TRIMDIR}"/*_1_trimmed.fastq.gz; do

  # Derive R2 and sample name from R1
  R2="${R1/_1_trimmed.fastq.gz/_2_trimmed.fastq.gz}"
  SAMPLE=$(basename "$R1" _1_trimmed.fastq.gz)
  OUTPATH="${OUTDIR}/${SAMPLE}"

  # Skip if already quantified
  if [[ -f "${OUTPATH}/quant.sf" ]]; then
    echo "[SKIP] ${SAMPLE} already quantified"
    continue
  fi

  # Verify R2 exists
  if [[ ! -f "$R2" ]]; then
    echo "ERROR: R2 file not found for ${SAMPLE}: expected ${R2}"
    exit 1
  fi

  echo "[$(date '+%H:%M:%S')] Quantifying ${SAMPLE}..."

  salmon quant \
    --index "$INDEX" \
    --libType A \
    -1 "$R1" \
    -2 "$R2" \
    --output "$OUTPATH" \
    --threads "$THREADS" \
    --gcBias \
    --validateMappings \
    --numBootstraps 30 \
    --seqBias \
    2> "${LOGDIR}/${SAMPLE}_salmon.log"

  # Extract and print mapping rate from log
  MAPRATE=$(grep "Mapping rate" "${LOGDIR}/${SAMPLE}_salmon.log" | tail -1 | awk '{print $NF}')
  echo "[$(date '+%H:%M:%S')] ✓ ${SAMPLE} — Mapping rate: ${MAPRATE}"

done

echo "[$(date '+%H:%M:%S')] All samples quantified."
```

Run it:

```bash
chmod +x scripts/bash/05_salmon_quant.sh
bash scripts/bash/05_salmon_quant.sh 2>&1 | tee logs/05_salmon_main.log
```

### What each flag does

| Flag | What it does | Why it matters |
|------|-------------|---------------|
| `--index` | Path to the Salmon index built in Part 1 | Required |
| `--libType A` | Auto-detect library strandedness | Safer than hard-coding; provides verification |
| `-1` / `-2` | Paired-end read files | Required for PE data — use `-r` for single-end |
| `--gcBias` | Correct for GC content bias in quantification | Important for CRC data — tumour samples often have GC biases |
| `--validateMappings` | More rigorous quasi-mapping that checks for false hits | Improves accuracy, small speed cost |
| `--numBootstraps 30` | Generate 30 bootstrap replicates of the counts | Enables downstream uncertainty quantification in `tximeta` |
| `--seqBias` | Correct for sequence-specific bias at read starts | TruSeq libraries have known hexamer priming bias |

> ⚠️ **`--numBootstraps 30`** increases runtime by approximately 3x. On a laptop this means ~10–15 minutes per sample instead of ~3 minutes. For 18 samples expect 3–4 hours total. Run this overnight or over a long break. The bootstraps are not strictly required for a standard DESeq2 analysis but they enable `tximeta` to propagate quantification uncertainty — a more rigorous approach.
>
> If you are short on time, remove `--numBootstraps 30` for now. You can always re-run later.

---

## 7. Checking mapping rates

The mapping rate is the single most important Salmon QC metric. It tells you what percentage of your trimmed reads were successfully mapped to the transcriptome.

**Expected for GSE50760:** 70–85% mapping rate. This dataset was originally aligned to hg19 with Tophat2; we are now aligning to GRCh38 with a newer transcriptome annotation, so rates may be slightly different from what the original authors reported.

| Mapping rate | Interpretation | Action |
|-------------|----------------|--------|
| >75% | Good | Proceed |
| 60–75% | Acceptable | Note in methods, investigate if consistently low |
| 50–60% | Concerning | Check species, check if trimming was too aggressive |
| <50% | Problematic | Do not proceed — investigate sample |

Run this after the main script completes to get a clean summary table:

```bash
# Print mapping rate summary for all samples
echo -e "Sample\tLibType\tMappingRate"
for LOG in logs/salmon/*_salmon.log; do
  SAMPLE=$(basename "$LOG" _salmon.log)
  LIBTYPE=$(grep "Automatically" "$LOG" | awk '{print $NF}' | tr -d '.')
  MAPRATE=$(grep "Mapping rate" "$LOG" | tail -1 | awk '{print $NF}')
  echo -e "${SAMPLE}\t${LIBTYPE}\t${MAPRATE}"
done
```

**Expected output:**

```
Sample          LibType   MappingRate
SRR1049158      IU        78.4321%
SRR1049159      IU        81.2045%
SRR1049160      IU        79.9832%
...
```

Two things to check in this output:

1. **All samples report `IU`** — if any sample reports a different library type that is a serious flag. It means either that sample had a different prep or there is a file swapping error.
2. **No sample below 60%** — investigate any outlier before proceeding to DESeq2.

> 💡 **Portfolio tip:** Include this mapping rate summary table in your GitHub README or final report. It demonstrates sample-level QC at the quantification stage, not just at the raw FASTQ stage.

---

## 8. Aggregating Salmon QC with MultiQC

MultiQC can read Salmon log files and produce a summary of mapping rates and library type detection across all samples — exactly like it did for FastQC in Part 2.

```bash
# Run MultiQC on Salmon output directories
multiqc results/salmon_quant/ \
  --outdir results/multiqc/salmon \
  --filename "GSE50760_salmon_multiqc" \
  --title "GSE50760 Salmon Quantification QC" \
  2>&1 | tee logs/06_multiqc_salmon.log
```

Open `results/multiqc/salmon/GSE50760_salmon_multiqc.html` in your browser. The key panels to look at:

**Salmon: Mapping rates** — bar chart of mapping rate per sample. Should be consistent across all 18 samples. Any sample more than 10 percentage points below the group median warrants investigation.

**Salmon: Observed library types** — should show all samples as `IU`. If any sample shows a different type it is highlighted automatically.

**Salmon: Fragment length distribution** — the distribution of inferred fragment lengths (insert sizes between paired reads). Should be a unimodal distribution peaking around 150–300 bp for TruSeq libraries. Bimodal distributions or very short peaks can indicate library quality issues.

---

## 9. Verifying output structure

After Salmon completes, your output directory should look like this:

```
results/salmon_quant/
├── SRR1049158/               ← normal colon, patient AMC_2
│   ├── quant.sf              ← transcript counts (this is what you need)
│   ├── lib_format_counts.json ← library type evidence
│   ├── aux_info/
│   │   ├── meta_info.json    ← run statistics including mapping rate
│   │   └── bootstrap/        ← bootstrap replicate count files
│   └── logs/
│       └── salmon_quant.log
├── SRR1049159/               ← primary CRC, patient AMC_2
│   └── quant.sf
├── SRR1049160/               ← liver metastasis, patient AMC_2
│   └── quant.sf
...
├── SRR1049175/               ← normal colon, patient AMC_24
├── SRR1049176/               ← primary CRC, patient AMC_24
└── SRR1049177/               ← liver metastasis, patient AMC_24
```

Run this verification to confirm all 18 samples completed successfully:

```bash
# Count quant.sf files — should be 18
find results/salmon_quant/ -name "quant.sf" | wc -l
# Expected: 18

# Check none are empty
for QSF in results/salmon_quant/*/quant.sf; do
  LINES=$(wc -l < "$QSF")
  if [[ "$LINES" -lt 1000 ]]; then
    echo "WARNING: ${QSF} has only ${LINES} lines — may be truncated"
  fi
done

# Preview one quant.sf file
head -5 results/salmon_quant/SRR1049158/quant.sf
```

**Expected output from head:**

```
Name                    Length  EffectiveLength  TPM         NumReads
ENST00000456328.2       1657    1508.000         0.000000    0.000
ENST00000450305.2       632     483.000          0.000000    0.000
ENST00000488147.1       1351    1202.000         4.823190    347.219
ENST00000619216.1       68      19.000           0.000000    0.000
```

Also verify the meta_info.json for one sample — this gives you the exact mapping rate as a number:

```bash
cat results/salmon_quant/SRR1049158/aux_info/meta_info.json | python3 -m json.tool | grep -E "percent_mapped|num_processed|num_mapped"
```

**Expected output:**

```json
"num_processed": 31204819,
"num_mapped": 24619441,
"percent_mapped": 78.897
```

---

## 10. Introduction to tximeta — what comes next

You now have 18 `quant.sf` files — one per sample. In Part 4 you will load these into R and build the DESeq2 input object. It is worth understanding the import step conceptually before you do it.

### Why tximeta, not direct import

You could read each `quant.sf` file into R manually with `read.table()` and assemble a count matrix yourself. This works but loses important information. `tximeta` is better for three reasons:

**1. Transcript-to-gene aggregation.** Salmon quantifies at transcript level — one row per transcript, ~230,000 rows. DESeq2 works at gene level — one row per gene, ~60,000 rows. `tximeta` collapses transcript counts to gene counts using the Ensembl annotation automatically, handling multi-transcript genes correctly.

**2. Provenance tracking.** `tximeta` records exactly which reference transcriptome was used, the Ensembl version, and the genome build — embedding this in the resulting R object. This makes your analysis reproducible and self-documenting.

**3. Uncertainty propagation.** If you ran Salmon with `--numBootstraps`, `tximeta` can import the bootstrap variance estimates alongside the point estimates. This information flows into downstream uncertainty-aware methods.

### What tximeta needs from you

Three things:

```r
# 1. A data frame linking sample names to quant.sf file paths
files <- file.path("results/salmon_quant", samples$sample_id, "quant.sf")

# 2. The sample metadata table (from data/metadata/sample_metadata.csv)
coldata <- read.csv("data/metadata/sample_metadata.csv")

# 3. The name of the Ensembl release used to build the Salmon index
# (Ensembl release 111 in our case — recorded in Part 1)
```

In Part 4 you will load all of this and build a `SummarizedExperiment` object — the R data structure that feeds directly into DESeq2. All the hard computational work is done. The rest is statistics and biology.

### The three-condition design

One thing to think about before Part 4: GSE50760 has **three conditions** — normal, primary, metastasis — from the **same 18 patients**. This is a matched triplet design. In DESeq2 this means:

- The design formula will be `~ patient_id + condition`
- You will run **three pairwise contrasts:**
  - Primary vs Normal — genes changing in early CRC
  - Metastasis vs Normal — genes changing across the full progression
  - Metastasis vs Primary — genes specifically associated with metastatic transition
- The metastasis vs primary contrast is the most scientifically interesting and clinically relevant

Thinking about this now means you arrive at Part 4 knowing what questions you are asking, not just running code.

---

## 11. Part 3 checkpoint

Before moving to Part 4, verify **all** of the following:

- [ ] Salmon ran successfully on all 18 samples
- [ ] `find results/salmon_quant/ -name "quant.sf" | wc -l` returns `18`
- [ ] All 18 samples show mapping rate >60% (check `logs/salmon/` or MultiQC report)
- [ ] All 18 samples report library type `IU` — auto-detected correctly
- [ ] No `quant.sf` file has fewer than 1,000 lines
- [ ] MultiQC Salmon report reviewed — no outlier samples
- [ ] Fragment length distribution is unimodal and peaks between 150–350 bp
- [ ] `meta_info.json` inspected for at least one sample to verify mapping rate numbers
- [ ] All changes committed to git

```bash
# Commit Part 3
git add scripts/bash/05_salmon_quant.sh results/salmon_quant/ logs/salmon/ results/multiqc/salmon/
git commit -m "Part 3 complete: Salmon quantification done, all 18 samples mapped"
```

---

## Coming up in Part 4

Part 4 is where the biology begins. You will move from the terminal into **R**, import all 18 `quant.sf` files using `tximeta`, build a DESeq2 dataset, run the three-condition paired design model, and produce your first differential expression results.

The key outputs of Part 4 will be:
- PCA plot showing separation of the three conditions
- Volcano plots for each of the three pairwise contrasts
- Ranked DEG tables saved as CSV files
- Identification of the top genes driving normal → primary → metastasis progression

| Part | Title | What you'll do |
|------|-------|----------------|
| ✅ Part 1 | Setup & Download | WSL, conda, SRA download, Salmon index |
| ✅ Part 2 | QC & Trimming | FastQC, MultiQC, Trimmomatic |
| ✅ Part 3 | Salmon Quantification | Pseudoalignment, quant.sf files, mapping QC |
| **▶ Part 4** | **DESeq2 Differential Expression** | **tximeta import, paired design, three contrasts, volcano plots** |
| Part 5 | GSEA & Pathway Analysis | Hallmark enrichment, metastasis gene signatures, figures |

---

*GSE50760 CRC RNA-seq Pipeline · Part 3 of 5 · Masters in Bioinformatics Portfolio Project*
