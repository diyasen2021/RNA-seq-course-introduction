# RNA-seq Pipeline
## Part 2 — Quality Control & Trimming

---

## Contents

1. [Why QC matters](#1-why-qc-matters)
2. [What FastQC measures](#2-what-fastqc-measures)
3. [Running FastQC on all 12 samples](#3-running-fastqc-on-all-12-samples)
4. [Interpreting FastQC reports](#4-interpreting-fastqc-reports)
5. [Aggregating with MultiQC](#5-aggregating-with-multiqc)
6. [What to look for across samples](#6-what-to-look-for-across-samples)
7. [Adapter trimming with Trimmomatic](#7-adapter-trimming-with-trimmomatic)
8. [Post-trimming QC](#8-post-trimming-qc)
9. [When to exclude a sample](#9-when-to-exclude-a-sample)
10. [Part 2 checkpoint](#10-part-2-checkpoint)

---

## 1. Why QC matters

The most important habit in bioinformatics is this: **never trust data you haven't looked at**.

Raw FASTQ files coming off a sequencer can have a range of quality problems — low-quality base calls at read ends, adapter sequences that weren't properly removed during library prep, unexpected GC content, PCR duplicates from over-amplified libraries, or RNA degradation artifacts. None of these problems announce themselves. If you skip QC and feed bad data straight into Salmon, the quantification will run without errors and produce plausible-looking numbers that are quietly wrong.

In cancer RNA-seq specifically, QC failures matter even more because:
- A single contaminated sample in a tumor/normal comparison can create hundreds of false positives in your DEG list
- Batch effects (different RNA extraction dates, different library prep kits) show up clearly in QC metrics before they destroy your PCA
- Poor-quality samples from certain TCGA collection sites are well-documented — you need to identify and handle them explicitly

QC is not a formality. It is where you build your first understanding of what your data actually looks like.

---

## 2. What FastQC measures

### What Does a FASTQ File Look Like?

A FASTQ file stores sequencing reads and their quality scores. Each read consists of **four lines**.

Example:

```text
@SRR123456.1
ATCGGATCGATCGATCGATC
+
IIIIIIIIIIIIIIIIIIIII
```

### Line 1: Sequence Identifier

```text
@SRR123456.1
```

This is the read header and contains information about the sequencing run, instrument, flowcell, lane, and read number.

### Line 2: Sequence

```text
ATCGGATCGATCGATCGATC
```

This is the nucleotide sequence determined by the sequencer.

### Line 3: Separator

```text
+
```

This line separates the sequence from the quality scores.

### Line 4: Quality Scores

```text
IIIIIIIIIIIIIIIIIIIII
```

Each character represents the quality score of the corresponding base in Line 2.

### The Phred quality score

The quality score Q is defined as:

```
Q = -10 × log₁₀(P_error)
```

| Q score | Error probability | Accuracy |
|---------|------------------|----------|
| Q10 | 1 in 10 | 90% |
| Q20 | 1 in 100 | 99% |
| Q30 | 1 in 1,000 | 99.9% |
| Q40 | 1 in 10,000 | 99.99% |

**Rule of thumb:** You want median per-base quality above Q28 across the entire read. Q20 is the minimum acceptable threshold — bases below Q20 should be trimmed.

# The Problem: FASTQ Stores Characters, Not Numbers

FASTQ files do not store quality scores as integers.

Instead, quality scores are encoded as ASCII characters.

Example:

```text
IIIIIIIIIIII
```

Software must decode these characters back into numerical quality scores.

---

# Phred+33 Encoding

Modern sequencing data uses:

```text
Quality Score + 33
```

| Q Score | ASCII Character |
|----------|----------------|
| 0 | ! |
| 10 | + |
| 20 | 5 |
| 30 | ? |
| 40 | I |

Example:

ASCII value of I = 73

Quality score:

```text
73 - 33 = 40
```

Therefore:

```text
I = Q40
```

---

# Phred+64 Encoding

Older Illumina instruments used:

```text
Quality Score + 64
```

| Q Score | ASCII Character |
|----------|----------------|
| 0 | @ |
| 10 | J |
| 20 | T |
| 30 | ^ |

Example:

ASCII value of T = 84

Quality score:

```text
84 - 64 = 20
```

Therefore:

```text
T = Q20
```

---

# Why Did Illumina Switch to Phred+33?

The old Phred+64 encoding caused several issues:

- Different sequencing platforms used different encodings.
- Bioinformatics tools frequently interpreted quality scores incorrectly.
- Data sharing between platforms was complicated.

Around 2011–2012, Illumina adopted **Phred+33**, which matched the original Sanger FASTQ specification.

---

# How FastQC Determines the Encoding

FastQC examines the range of quality score characters present in the FASTQ file.

If it sees characters such as:

```text
!
"
#
$
%
```

the file must be Phred+33 because these characters cannot occur in Phred+64 encoding.

Typical FastQC output:

```text
Encoding: Sanger / Illumina 1.9
```

or

```text
Encoding: Illumina 1.9
```

Both indicate Phred+33 encoding.

---

# Modern RNA-seq and FASTQ Files

For data generated on modern Illumina platforms such as:

- NovaSeq
- NextSeq
- MiSeq

you can safely assume the FASTQ files use **Phred+33** encoding.

Most contemporary bioinformatics tools expect Phred+33 by default.

---

### What does FastQC measure?
FastQC reads your FASTQ file and runs 11 diagnostic modules. Each module either **passes** (green ✅), **warns** (orange ⚠️), or **fails** (red ❌). Importantly, these flags are calibrated for generic data — some "failures" are expected and acceptable in RNA-seq. You need to interpret them in context, not just count green ticks.

| Module | What it measures | Expected in RNA-seq |
|--------|-----------------|---------------------|
| **Per-base sequence quality** | Phred quality score at each position across all reads | High (Q30+) at 5' end, slight drop at 3' end is normal |
| **Per-tile sequence quality** | Quality by flow cell tile — detects localised sequencer problems | Should be uniform; bright red tiles = sequencer issue |
| **Per-sequence quality scores** | Distribution of mean quality per read | Sharp peak at Q35–40 is ideal |
| **Per-base sequence content** | % A/T/G/C at each position | First 10–15 bases often biased in RNA-seq — normal |
| **Per-sequence GC content** | GC distribution across reads vs. theoretical | Slight deviation is normal; sharp peaks = contamination |
| **Per-base N content** | Frequency of uncalled bases (N) | Should be near zero throughout |
| **Sequence length distribution** | Read length distribution | Should be a single sharp peak at your read length |
| **Sequence duplication levels** | % of reads that are exact duplicates | RNA-seq: 20–50% duplication is acceptable; >80% is concerning |
| **Overrepresented sequences** | Sequences appearing far more than expected | Adapters, rRNA, poly-A tails — identifies contamination |
| **Adapter content** | Illumina adapter sequences found in reads | Should be near zero after good library prep; if high, trimming is needed |
| **Kmer content** | Enriched k-mers at specific positions | Often fails in RNA-seq — usually not actionable |



## 3. Running FastQC on all 12 samples

Save this as `scripts/bash/02_fastqc_pretrim.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────
INDIR="data/raw_fastq"
OUTDIR="results/fastqc/pre_trim"
THREADS=6
mkdir -p "$OUTDIR"

echo "[$(date '+%H:%M:%S')] Running FastQC on all samples..."

# Run FastQC on all FASTQ files in parallel
fastqc \
  "${INDIR}"/*.fastq.gz \
  --outdir "$OUTDIR" \
  --threads "$THREADS" \
  --extract \          # Unzip the output HTML — makes it easier to browse
  2>&1 | tee logs/02_fastqc_pretrim.log

echo "[$(date '+%H:%M:%S')] FastQC complete. Reports in ${OUTDIR}"
```

Run it:

```bash
chmod +x scripts/bash/02_fastqc_pretrim.sh
bash scripts/bash/02_fastqc_pretrim.sh
```

**Expected output:**
```
[09:14:22] Running FastQC on all samples...
Started analysis of SRR1616919.fastq.gz
Approx 5% complete for SRR1616919.fastq.gz
...
Analysis complete for SRR1616930.fastq.gz
[09:31:05] FastQC complete. Reports in results/fastqc/pre_trim
```

FastQC generates two files per sample:
- `SRR1616919_fastqc.html` — open this in a browser to view the full report
- `SRR1616919_fastqc.zip` — contains the raw data behind the plots (used by MultiQC)

> ⏱️ **Runtime:** FastQC takes approximately 2–5 minutes per sample. With `--threads 6` and 12 samples, expect 5–15 minutes total.

---

## 4. Interpreting FastQC reports

Open any `_fastqc.html` file in your browser. Here is what to look for in each module for TCGA-KIRC data specifically.

### 4.1 Per-base sequence quality

This is the most important plot. It shows a box-and-whisker plot of Phred quality at each base position.

**What good looks like:**
```
Quality
40 |████████████████████████████████████▓▓▓▓▓▓▓░░
30 |─────────────────────────────────────────────
20 |
    1   10   20   30   40   50  (base position)
```
- Median (orange line) stays above Q28 for the full read length
- Lower whiskers stay above Q20
- A gradual decline at the 3' end (last 10–15 bases) is normal — sequencing chemistry degrades slightly

**Warning signs:**
- Median drops below Q28 before position 40 → aggressive trimming needed
- Sudden drop at a specific position in all samples → sequencer issue
- Very low quality throughout → library prep problem, sample may need exclusion

### 4.2 Per-base sequence content

RNA-seq almost always "fails" this module and that is **expected and acceptable**. The first 10–15 bases of RNA-seq reads show biased nucleotide composition because of:
- Random hexamer priming bias during cDNA synthesis
- Sequence-specific priming preferences

**What to ignore:** Biased composition in the first 10–15 bases.
**What to flag:** Biased composition in the middle or end of reads, which can indicate contamination.

### 4.3 Sequence duplication levels

This module flags frequently if duplication exceeds 20%, which will happen with TCGA data. This is **not a problem** for RNA-seq:

- Highly expressed genes (ACTB, GAPDH, ribosomal genes) produce enormous numbers of nearly identical reads
- Unlike DNA-seq where duplication indicates PCR over-amplification, in RNA-seq it reflects genuine biological abundance
- **Do not deduplicate RNA-seq data** — you will throw away real signal

> ⚠️ A duplication warning in RNA-seq does not mean your data is bad. It means some genes are very highly expressed. This is expected.

### 4.4 Adapter content

This is the module that tells you whether trimming is necessary. Illumina TruSeq adapters begin appearing in reads when the cDNA insert is shorter than the read length — the sequencer reads through the insert and into the adapter sequence.

**What you are looking for:**
```
% adapter
15 |          ░░░░░░░░░░░░░░░░░▓▓▓▓▓▓▓████████
10 |      ░░░
 5 |
 0 |████
    1   10   20   30   40   50  (base position)
```
- Near zero throughout → library is clean, trimming may still be good practice
- Rising curve from position 30+ → adapter contamination present, trimming is required
- High adapter content from position 1 → very short inserts, check library prep

### 4.5 Overrepresented sequences

FastQC will flag sequences appearing in >0.1% of reads. For TCGA RNA-seq, common culprits are:

| Sequence type | What it means | Action |
|---------------|---------------|--------|
| Illumina TruSeq adapter | Read-through into adapter | Trim with Trimmomatic |
| Poly-A tail | Normal RNA-seq artefact from mRNA poly-A selection | Acceptable, soft-trim if excessive |
| rRNA sequences | Incomplete ribo-depletion | Flag — check if consistent across samples |
| No overrepresented sequences | Clean library | ✅ Ideal |

---

## 5. Aggregating with MultiQC

Looking at 12 individual FastQC HTML files is tedious and makes it hard to spot patterns across samples. MultiQC reads all the FastQC output files and produces a single interactive HTML report.

```bash
# Run MultiQC on the FastQC output directory
multiqc results/fastqc/pre_trim/ \
  --outdir results/multiqc/pre_trim \
  --filename "KIRC_pre_trim_multiqc" \
  --title "TCGA-KIRC Pre-Trimming QC" \
  2>&1 | tee logs/03_multiqc_pretrim.log
```

Open `results/multiqc/pre_trim/KIRC_pre_trim_multiqc.html` in your browser.

**Expected output:**
```
[INFO   ]         multiqc : This is MultiQC v1.14
[INFO   ]         multiqc : Template    : default
[INFO   ]          fastqc : Found 12 reports
[INFO   ]         multiqc : Report      : results/multiqc/pre_trim/KIRC_pre_trim_multiqc.html
[INFO   ]         multiqc : Data        : results/multiqc/pre_trim/KIRC_pre_trim_multiqc_data
```

---

## 6. What to look for across samples

The MultiQC report lets you compare all 12 samples side by side. Here is what to examine systematically.

### 6.1 General statistics table

The first table in MultiQC shows key metrics per sample. Check:

| Metric | Acceptable range | Flag if |
|--------|-----------------|---------|
| Total reads | 20M–50M | <15M per sample |
| % GC | 45–60% | >65% or <40% |
| % Dups | 20–70% | >85% (RNA-seq) |
| % Failed | <5% | >20% reads failed QC |
| Median quality | >Q28 | <Q25 |

### 6.2 Sequence quality heatmap

MultiQC plots mean quality per position as a heatmap across all samples. You are looking for:
- **Consistent pattern across samples** → sequencing run was uniform
- **One or two outlier samples with much lower quality** → those samples are candidates for exclusion
- **All samples dropping at the same position** → systematic issue with this sequencing run

### 6.3 Adapter content panel

This panel overlays adapter content curves for all 12 samples. In TCGA data you may see:
- Low-level adapter content (1–5%) → trim as standard practice
- Higher adapter content (5–15%) in some samples → those samples had shorter inserts

### 6.4 Per-sample status summary

At the bottom of the MultiQC report is a pass/warn/fail summary grid for every module. Use this to quickly identify which samples need closer inspection.

> 💡 **Portfolio tip:** Take a screenshot of the MultiQC summary table and include it in your GitHub README or final report. It demonstrates that you performed QC rigorously — not just that you ran the tool.

---

## 7. Adapter trimming with Trimmomatic

Even if adapter content looks low, trimming is standard practice and costs nothing. Trimmomatic removes adapter sequences and low-quality bases from the 3' end of reads.

### 7.1 Understanding Trimmomatic parameters

Trimmomatic processes reads through a series of steps in the order you specify them:

| Step | Parameter | What it does |
|------|-----------|-------------|
| Adapter removal | `ILLUMINACLIP` | Removes Illumina adapter sequences |
| Leading quality | `LEADING:3` | Removes bases from 5' end with quality < 3 |
| Trailing quality | `TRAILING:20` | Removes bases from 3' end with quality < 20 |
| Sliding window | `SLIDINGWINDOW:4:20` | Cuts when mean quality in a 4-base window drops below Q20 |
| Minimum length | `MINLEN:36` | Discards reads shorter than 36 bases after trimming |

### 7.2 Adapter file

Trimmomatic needs to know which adapter sequences to look for. The standard Illumina TruSeq adapters are bundled with Trimmomatic:

```bash
# Find the adapter file in your conda environment
find $CONDA_PREFIX -name "TruSeq3-SE.fa" 2>/dev/null
# Should return something like:
# /home/user/miniforge3/envs/rnaseq/share/trimmomatic/adapters/TruSeq3-SE.fa
```

> **SE vs PE adapters:** TCGA-KIRC SRA data is single-end, so use `TruSeq3-SE.fa`. If you ever work with paired-end data use `TruSeq3-PE.fa`. Using the wrong adapter file will result in no adapter removal — a silent failure.

### 7.3 Trimming script

Save this as `scripts/bash/03_trimmomatic.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────
INDIR="data/raw_fastq"
OUTDIR="data/trimmed"
LOGDIR="logs/trimmomatic"
THREADS=6
mkdir -p "$OUTDIR" "$LOGDIR"

# Locate adapter file bundled with Trimmomatic in conda environment
ADAPTERS=$(find "$CONDA_PREFIX" -name "TruSeq3-SE.fa" | head -1)

if [[ -z "$ADAPTERS" ]]; then
  echo "ERROR: Could not find TruSeq3-SE.fa adapter file"
  echo "Check your conda environment: find \$CONDA_PREFIX -name '*.fa' | grep -i adapter"
  exit 1
fi

echo "[$(date '+%H:%M:%S')] Using adapters: $ADAPTERS"

# ── Trimming loop ───────────────────────────────────────────────────────
for FASTQ in "${INDIR}"/*.fastq.gz; do

  SAMPLE=$(basename "$FASTQ" .fastq.gz)
  OUTFILE="${OUTDIR}/${SAMPLE}_trimmed.fastq.gz"

  # Skip if already trimmed
  if [[ -f "$OUTFILE" ]]; then
    echo "[SKIP] ${SAMPLE} already trimmed"
    continue
  fi

  echo "[$(date '+%H:%M:%S')] Trimming ${SAMPLE}..."

  trimmomatic SE \
    -threads "$THREADS" \
    -phred33 \
    "$FASTQ" \
    "$OUTFILE" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
    LEADING:3 \
    TRAILING:20 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    2> "${LOGDIR}/${SAMPLE}_trimmomatic.log"

  # Print summary from log
  echo "[$(date '+%H:%M:%S')] ✓ ${SAMPLE} — $(grep 'Reads Written' ${LOGDIR}/${SAMPLE}_trimmomatic.log)"

done

echo "[$(date '+%H:%M:%S')] All samples trimmed."
```

Run it:

```bash
chmod +x scripts/bash/03_trimmomatic.sh
bash scripts/bash/03_trimmomatic.sh 2>&1 | tee logs/03_trimmomatic_main.log
```

**Expected output per sample:**
```
[10:02:14] Trimming SRR1616919...
[10:04:31] ✓ SRR1616919 — Reads Written (surviving): 27,844,103 (97.94%)
[10:04:31] Trimming SRR1616920...
[10:06:58] ✓ SRR1616920 — Reads Written (surviving): 30,621,447 (98.13%)
```

> ⚠️ **Flag for exclusion:** If a sample shows fewer than 90% reads surviving, investigate that sample's pre-trim FastQC report before proceeding. A sample losing >15% of reads to trimming usually has a library quality problem.

### 7.4 Understanding the Trimmomatic log

Each sample produces a log file in `logs/trimmomatic/`. It contains a summary like:

```
Input Reads: 28431207
Surviving: 27844103 (97.94%)
Dropped: 587104 (2.06%)
```

Record these numbers — you will reference them in your methods section and they are part of your QC documentation for the portfolio project.

---

## 8. Post-trimming QC

Always run FastQC again after trimming to confirm that:
1. Adapter content is now gone (or substantially reduced)
2. Per-base quality has improved at read ends
3. Read length distribution has shifted slightly shorter (expected — trimming removes bases)
4. You have not accidentally over-trimmed (losing too many reads or making reads very short)

```bash
#!/usr/bin/env bash
# scripts/bash/04_fastqc_posttrim.sh
set -euo pipefail

INDIR="data/trimmed"
OUTDIR="results/fastqc/post_trim"
THREADS=6
mkdir -p "$OUTDIR"

echo "[$(date '+%H:%M:%S')] Running post-trimming FastQC..."

fastqc \
  "${INDIR}"/*_trimmed.fastq.gz \
  --outdir "$OUTDIR" \
  --threads "$THREADS" \
  --extract \
  2>&1 | tee logs/04_fastqc_posttrim.log

# Run MultiQC on post-trim reports
multiqc "$OUTDIR" \
  --outdir results/multiqc/post_trim \
  --filename "KIRC_post_trim_multiqc" \
  --title "TCGA-KIRC Post-Trimming QC" \
  2>&1 | tee logs/04_multiqc_posttrim.log

echo "[$(date '+%H:%M:%S')] Post-trim QC complete."
```

```bash
bash scripts/bash/04_fastqc_posttrim.sh
```

### Comparing pre- and post-trimming

Open both MultiQC reports side by side and compare:

| Module | Pre-trim | Post-trim (expected) |
|--------|----------|---------------------|
| Per-base quality | May dip at 3' end | Improved at 3' end |
| Adapter content | May show adapter signal | Near zero throughout |
| Sequence length | Sharp peak at 50bp | Broader distribution, peak slightly lower |
| Total reads | Baseline | 95–99% of pre-trim (good trimming loses very few reads) |

> 💡 **Portfolio tip:** Include both pre- and post-trim MultiQC screenshots in your report or README. The before/after comparison demonstrates that your trimming was effective and conservative — not aggressive. A reviewer can immediately see that your trimming improved quality without destroying reads.

---

## 9. When to exclude a sample

Most of your 12 samples will pass QC without issues. But you need to know the criteria for exclusion in case they don't. A sample should be excluded if it shows **two or more** of the following:

| Problem | Threshold | Module |
|---------|-----------|--------|
| Very low read count | <15M reads after trimming | General stats |
| Poor overall quality | Median per-base Q < 25 throughout | Per-base quality |
| Extreme GC content | <40% or >65% | Per-sequence GC |
| Massive adapter contamination | >20% adapter at position 30 | Adapter content |
| Very high duplication | >85% (in the context of low read count) | Duplication |
| rRNA contamination | Overrepresented sequences match rRNA | Overrepresented sequences |

> **Important:** Excluding a sample requires removing it from your `sample_metadata.csv` as well. Downstream tools (DESeq2) will fail or produce incorrect results if your metadata and count files are mismatched.

If you do need to exclude a sample:

```bash
# Document the exclusion in your log
echo "EXCLUDED: SRR1616XXX — reason: <describe problem>" >> logs/sample_exclusions.txt

# Remove from metadata
# Edit data/metadata/sample_metadata.csv to remove that row

# Commit the change with a clear message
git add data/metadata/sample_metadata.csv logs/sample_exclusions.txt
git commit -m "Exclude SRR1616XXX: failed QC due to <reason>"
```

Documenting exclusions explicitly — rather than silently removing samples — is a mark of rigorous, reproducible analysis.

---
