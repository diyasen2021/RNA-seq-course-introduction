RNA-seq Pipeline
## Part 1 — Environment Setup, WSL Installation & Data Download

---

## Contents

1. [Scientific background & project overview](#1-scientific-background--project-overview)
2. [System requirements & software stack](#2-system-requirements--software-stack)
3. [Installing WSL2 (Windows users)](#3-installing-wsl2-windows-users)
4. [Installing bioinformatics tools](#4-installing-bioinformatics-tools)
5. [Setting up your project directory](#5-setting-up-your-project-directory)
6. [Downloading TCGA-KIRC metadata & SRA accessions](#6-downloading-tcga-kirc-metadata--sra-accessions)
7. [Downloading FASTQ files with fasterq-dump](#7-downloading-fastq-files-with-fasterq-dump)
8. [Verifying your downloads](#8-verifying-your-downloads)
9. [Part 1 checkpoint](#9-part-1-checkpoint)

---

## 1. Scientific background & project overview

Colorectal cancer (CRC) is the third most common cancer worldwide and the second leading cause of cancer-related death. While early-stage CRC is highly treatable, the development of liver metastasis — which occurs in approximately 50% of patients — dramatically worsens prognosis, with 5-year survival dropping from ~90% for localised disease to ~14% for metastatic disease.
The molecular mechanisms driving this metastatic transition are incompletely understood, making it one of the most active areas of cancer transcriptomics research. RNA-seq allows us to ask directly: what changes in gene expression enable a primary colorectal tumour to invade, survive in circulation, and colonise the liver?
This makes GSE50760 an exceptional teaching dataset:

The three-condition design (normal → primary tumour → liver metastasis) captures disease progression, not just presence or absence of cancer
All three samples come from the same 18 patients — a perfectly matched triplet design that controls for inter-patient genetic variation
The biology is clinically urgent — metastasis is what kills CRC patients
The dataset has been cited and reanalysed multiple times in the literature, so your results are directly benchmarkable

- The three-condition design (normal → primary tumour → liver metastasis) captures disease **progression**, not just presence or absence of cancer
- All three samples come from the **same 18 patients** — a perfectly matched triplet design that controls for inter-patient genetic variation
- The biology is clinically urgent — metastasis is what kills CRC patients
- The dataset has been cited and reanalysed multiple times in the literature, so your results are directly benchmarkable

### What you will produce by the end of this project

| Analysis | Description |
|----------|-------------|
| **Differential Expression (3 contrasts)** | Primary vs Normal, Metastasis vs Normal, and Metastasis vs Primary — identifying genes driving each stage of progression |
| **Metastasis gene signature** | Genes specifically upregulated in liver metastasis but not primary tumour — the transcriptional fingerprint of metastatic transition |
| **Pathway Enrichment (GSEA)** | GSEA against Hallmark and KEGG gene sets. Expected hits: `EPITHELIAL_MESENCHYMAL_TRANSITION`, `ANGIOGENESIS`, `TGF_BETA_SIGNALLING` |
| **Visualisation & reporting** | PCA, volcano plots, and heatmaps publication-ready enough for a portfolio or paper figure |

> **Note on dataset scope:** In this tutorial you will download and process a **subset of 18 samples** — one triplet (normal, primary, metastasis) from each of 6 patients — through the full FASTQ-to-counts pipeline to learn every step hands-on. In Part 4 (DESeq2), you will load the **full GSE50760 cohort** (all 54 samples, all 18 patients) using the NCBI-generated raw counts matrix directly from GEO. This mirrors real-world practice: learn the pipeline on a manageable subset, do the science on the complete dataset.

---

## 2. System requirements & software stack

### Hardware requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| RAM | 8 GB | 16 GB | Salmon needs ~8 GB for the human transcriptome index. STAR would need 32 GB — this is why we use Salmon. |
| Disk | 80 GB free | 150 GB free | 12 paired-end FASTQ files ≈ 60–80 GB. Intermediate files add ~20 GB. Full cohort counts matrix is <50 MB. |
| CPU | 4 cores | 8+ cores | Salmon is highly parallelized. More threads = faster quantification. |
| OS | Windows 10 (WSL2) | Ubuntu 22.04 LTS / macOS | All pipeline tools run natively on Linux. WSL2 on Windows works well. |
| Internet | Stable broadband | 50+ Mbps | Downloading 12 FASTQ files will transfer ~60–80 GB from NCBI SRA servers. |

### Software stack

| Tool | Version | Purpose | Required? |
|------|---------|---------|-----------|
| SRA Toolkit | 3.0+ | Download FASTQ files from NCBI SRA | ✅ Required |
| FastQC | 0.12+ | Per-sample quality control reports | ✅ Required |
| MultiQC | 1.14+ | Aggregate QC reports across all samples | ✅ Required |
| Trimmomatic | 0.39+ | Adapter removal and quality trimming | ✅ Required |
| Salmon | 1.10+ | Pseudoalignment + transcript quantification | ✅ Required |
| R | 4.3+ | All downstream statistical analysis | ✅ Required |
| Conda / Mamba | latest | Package management | ⭐ Recommended |

---

## 3. Installing WSL2 (Windows users)

> **macOS / Linux users:** Skip this section entirely. Open your terminal and jump to [Section 4](#4-installing-bioinformatics-tools).

WSL2 (Windows Subsystem for Linux 2) gives you a full Ubuntu environment inside Windows. All bioinformatics tools will be installed inside this Linux environment.

### 3.1 Enable WSL2

Open **PowerShell as Administrator** (right-click Start → Windows PowerShell (Admin)) and run:

```powershell
# Enable WSL and install Ubuntu 22.04 in one command (Windows 11 / Win10 build 19041+)
wsl --install

# If that doesn't work on older Windows 10, run these separately:
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```

> ⚠️ **Restart required.** You must restart your computer after running these commands before continuing.

### 3.2 Set WSL2 as default and install Ubuntu

```powershell
# After restart — set WSL version 2 as default
wsl --set-default-version 2

# Install Ubuntu 22.04 LTS (if not already installed)
wsl --install -d Ubuntu-22.04

# Verify WSL2 is active
wsl --list --verbose
```

**Expected output:**
```
  NAME            STATE           VERSION
* Ubuntu-22.04    Running         2
```

When Ubuntu launches for the first time it will ask you to create a Unix username and password. This is separate from your Windows login. From now on, all pipeline commands run inside the Ubuntu terminal.

### 3.3 Update the base system

Open the Ubuntu app (search "Ubuntu" in Start menu) and run:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y curl wget git build-essential unzip
```

### 3.4 Configure WSL2 memory (important for Salmon)

By default WSL2 may not use all your RAM. Create a config file to give it more.

Run in PowerShell (not Ubuntu):
```powershell
notepad "$env:USERPROFILE\.wslconfig"
```

In the Notepad window that opens, paste the following — adjust numbers to your actual RAM, leaving ~4 GB for Windows:

```ini
[wsl2]
memory=12GB        # Use 12GB if you have 16GB RAM; 6GB if you have 8GB
processors=6       # Use n-2 of your CPU cores
swap=8GB           # Extra virtual memory — useful for Salmon index build
localhostForwarding=true
```

Save, close Notepad, then restart WSL2 from PowerShell:

```powershell
wsl --shutdown
# Wait 5 seconds, then reopen the Ubuntu terminal
```

---

## 4. Installing bioinformatics tools

We will use **Miniforge** (a minimal conda installer) to manage all pipeline tools. This keeps everything isolated and reproducible — a critical habit for portfolio-quality work.

### 4.1 Install Miniforge

```bash
# Download Miniforge installer
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Run installer (answer yes to all prompts)
bash Miniforge3-Linux-x86_64.sh

# Reload shell to activate conda
source ~/.bashrc

# Verify installation
conda --version
# conda 23.11.0
```

### 4.2 Create the pipeline environment

Create a dedicated conda environment with all required tools. Using a single environment means the exact software versions are captured and reproducible — include this in your portfolio README.

```bash
# Create environment named 'rnaseq' with all pipeline tools
conda create -n rnaseq -c bioconda -c conda-forge \
    sra-tools fastqc multiqc trimmomatic salmon \
    r-base=4.3 bioconductor-deseq2 bioconductor-tximeta \
    r-tidyverse r-ggplot2 r-pheatmap -y

# Activate the environment (run this every time you open a new terminal)
conda activate rnaseq

# Verify key tools are installed
fasterq-dump --version
# fasterq-dump version 3.0.10

salmon --version
# salmon 1.10.3

fastqc --version
# FastQC v0.12.1
```

> 💡 **Portfolio tip — export your environment file.** After installation, run the following and commit the output to your GitHub repo. This lets anyone reproduce your exact software stack:
> ```bash
> conda env export > environment.yml
> git add environment.yml
> git commit -m "Add conda environment file"
> ```

---

## 5. Setting up your project directory

Organised directory structure is not optional for portfolio work — it signals professionalism and makes your project reproducible. We follow a convention used in production bioinformatics pipelines.

### 5.1 Create the directory tree

```bash
# Create the full project directory tree
mkdir -p ~/kirc_rnaseq/{
    data/raw_fastq,
    data/trimmed,
    data/metadata,
    reference/transcriptome,
    reference/salmon_index,
    results/fastqc,
    results/multiqc,
    results/salmon_quant,
    results/deseq2,
    results/gsea,
    results/figures,
    scripts/bash,
    scripts/R,
    logs
}

cd ~/kirc_rnaseq
```

Your project structure should look like this:

```
kirc_rnaseq/
├── data/
│   ├── raw_fastq/       ← FASTQ files downloaded from SRA
│   ├── trimmed/         ← Trimmomatic output (Part 2)
│   └── metadata/        ← SRA accession list, clinical metadata
├── reference/
│   ├── transcriptome/   ← Ensembl cDNA FASTA (downloaded below)
│   └── salmon_index/    ← Salmon index built from transcriptome
├── results/
│   ├── fastqc/          ← FastQC HTML reports per sample
│   ├── multiqc/         ← Aggregated MultiQC report
│   ├── salmon_quant/    ← Salmon quant.sf files per sample
│   ├── deseq2/          ← DESeq2 results tables & RDS objects
│   ├── gsea/            ← GSEA enrichment results
│   └── figures/         ← All publication-quality plots
├── scripts/
│   ├── bash/            ← Shell scripts: download, trim, quantify
│   └── R/               ← R scripts: DESeq2, GSEA, survival
├── logs/                ← stdout/stderr from each pipeline step
├── environment.yml      ← Conda environment (committed to git)
└── README.md            ← Your project documentation
```

### 5.2 Initialise a git repository

Start version-controlling your project now. Commit after each major step.

```bash
git init

# Create a .gitignore — do NOT commit raw FASTQ files (too large)
cat > .gitignore << 'EOF'
data/raw_fastq/
data/trimmed/
reference/transcriptome/
*.fastq.gz
*.bam
*.sra
.snakemake/
EOF

# Create README
cat > README.md << 'EOF'
# TCGA-KIRC RNA-seq Pipeline
Full RNA-seq pipeline: FASTQ → Salmon → DESeq2 → GSEA → Survival analysis
Dataset: TCGA Kidney Renal Clear Cell Carcinoma (KIRC)
Samples: 6 tumor + 6 matched normal (pipeline) + full cohort via recount3
EOF

git add .gitignore README.md
git commit -m "Initial project structure"
```

---

## 6.  Get the SRR accession list from SRA Run Selector

Go to the SRA Run Selector for this dataset: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP029880

On that page:

Click "Metadata" — downloads SraRunTable.txt, a CSV with one row per sample containing SRR accession, sample name, condition, patient ID, and sequencing metadata.
Click "Accession List" — downloads SRR_Acc_List.txt, a plain text file with one SRR accession per line.

Download both files and move them into your project:

```bash
mv ~/Downloads/SraRunTable.txt data/metadata/SraRunTable.txt
mv ~/Downloads/SRR_Acc_List.txt data/metadata/SRR_Acc_List.txt
```

### 6.2 Inspect the SraRunTable

This CSV file is what connects sample IDs to biological conditions in DESeq2. Check what columns are available


```bash
head -2 data/metadata/SraRunTable.txt | tr ',' '\n' | nl
```

Select your 18-sample subset (6 patients × 3 conditions). For the pipeline tutorial we want 6 patients, each contributing all three samples — normal, primary, metastasis. This preserves the matched triplet structure.Then create your subset accession list. Save it in SRR_accessions.txt.

## 7. Downloading FASTQ files with fasterq-dump

`fasterq-dump` is the modern, parallelised replacement for `fastq-dump`. It is significantly faster and handles paired-end splitting correctly by default. We write a bash script rather than running commands manually because scripted downloads are reproducible and resumable.

### 7.1 Configure SRA Toolkit

Before downloading, configure the SRA toolkit to use your project directory as cache:

```bash
# Configure non-interactively
vdb-config --set /repository/user/main/public/root=$HOME/kirc_rnaseq/data/raw_fastq

# Or run the interactive wizard
vdb-config --interactive
# Go to CACHE tab → set location to ~/kirc_rnaseq/data/raw_fastq/
# Press 's' to save, 'x' to exit
```

### 7.2 Download script

Save this as `scripts/bash/01_download_fastq.sh`. The script logs progress, handles errors, and can be restarted if interrupted — already-downloaded files are automatically skipped.

```bash
#!/usr/bin/env bash
set -euo pipefail  # Exit on error, undefined vars, pipe failures

# ── Configuration ──────────────────────────────────────────────────────
ACCESSIONS="data/metadata/SRR_accessions.txt"
OUTDIR="data/raw_fastq"
THREADS=6          # Adjust to your CPU core count
TMPDIR="/tmp/sra_tmp"
mkdir -p "$OUTDIR" "$TMPDIR"

# ── Download loop ───────────────────────────────────────────────────────
while read -r SRR TYPE PATIENT || [[ -n "$SRR" ]]; do

  # Skip comment lines and empty lines
  [[ "$SRR" =~ ^# ]] && continue
  [[ -z "$SRR" ]] && continue

  # Skip if already downloaded
  if [[ -f "${OUTDIR}/${SRR}.fastq.gz" ]]; then
    echo "[SKIP] ${SRR} already exists"
    continue
  fi

  echo "[$(date '+%H:%M:%S')] Downloading ${SRR} (${TYPE}, ${PATIENT})"

  fasterq-dump "$SRR" \
    --outdir "$OUTDIR" \
    --threads "$THREADS" \
    --temp "$TMPDIR" \
    --split-files \       # Creates _1.fastq and _2.fastq for paired-end
    --skip-technical \    # Exclude technical reads (barcodes etc.)
    --progress            # Show download progress bar

  # Compress immediately to save disk space
  echo "[$(date '+%H:%M:%S')] Compressing ${SRR}..."
  gzip "${OUTDIR}/${SRR}"*.fastq

  echo "[$(date '+%H:%M:%S')] ✓ ${SRR} complete"

done < "$ACCESSIONS"

echo "[$(date '+%H:%M:%S')] All downloads complete."
```

Run it:

```bash
chmod +x scripts/bash/01_download_fastq.sh
bash scripts/bash/01_download_fastq.sh 2>&1 | tee logs/01_download.log
```

### 7.3 Download the Salmon reference transcriptome

Salmon maps reads against the **transcriptome** (cDNA), not the full genome. This requires only ~8 GB RAM versus 32 GB for genome alignment.

```bash
# Download Ensembl GRCh38 cDNA reference
wget -P reference/transcriptome/ \
  https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Also download the GTF for gene-level annotation
wget -P reference/ \
  https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz

# Build the Salmon index (~10 minutes, uses ~8 GB RAM at peak)
salmon index \
  -t reference/transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz \
  -i reference/salmon_index \
  --gencode \
  -p 6 \
  2>&1 | tee logs/salmon_index.log
```

**Expected output (last few lines):**
```
[info] Building new index...
[info] done building the index
[info] Index building took 8 minutes 23 seconds
```

---

## 8. Verifying your downloads

Never proceed to QC or alignment without verifying your downloads. Truncated FASTQ files cause cryptic errors later in the pipeline.

```bash
# Check all 12 files exist and have non-zero size
ls -lh data/raw_fastq/*.fastq.gz

# Count reads in each file (a valid 30M-read sample has exactly 120M lines)
for f in data/raw_fastq/*.fastq.gz; do
  echo -n "$f: "
  zcat "$f" | wc -l | awk '{print $1/4, "reads"}'
done

# Verify gzip integrity — detects truncation
for f in data/raw_fastq/*.fastq.gz; do
  gzip -t "$f" && echo "OK: $f" || echo "CORRUPTED: $f"
done
```

**Expected output:**
```
data/raw_fastq/SRR1616919.fastq.gz: 28431207 reads
data/raw_fastq/SRR1616920.fastq.gz: 31204819 reads
data/raw_fastq/SRR1616921.fastq.gz: 25918430 reads
... (each file should show 20M–45M reads)

OK: data/raw_fastq/SRR1616919.fastq.gz
OK: data/raw_fastq/SRR1616920.fastq.gz
...
```

> ❌ **If a file shows "CORRUPTED" or fewer than 5M reads:** Delete the file and re-run the download script — it will skip intact files and re-download only the failed ones:
> ```bash
> rm data/raw_fastq/SRR1616XXX.fastq.gz
> bash scripts/bash/01_download_fastq.sh
> ```

### Commit your progress

```bash
# Commit scripts and metadata (not FASTQ files — they're in .gitignore)
git add scripts/ data/metadata/ logs/ reference/salmon_index/ environment.yml
git commit -m "Part 1 complete: data downloaded, Salmon index built"
```

---



*TCGA-KIRC RNA-seq Pipeline · Part 1 of 5 · Masters in Bioinformatics Portfolio Project*
