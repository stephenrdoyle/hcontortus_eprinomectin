# RNA-seq analysis pipeline

This repository describes the RNA-seq workflow used to process paired-end transcriptomic sequencing data from *Haemonchus contortus*, from raw-read acquisition through to the generation of a raw gene-count matrix.

Reads were aligned against the *H. contortus* PRJEB506 WBPS18 reference genome, and gene-level counts were obtained using the corresponding WBPS18 canonical gene annotation.

---

# Workflow overview

```text
Paired-end FASTQ
       │
       ▼
0. Read acquisition
       │
       ▼
1. Raw-read quality control
   FastQC
       │
       ▼
2. QC summary
   MultiQC
       │
       ▼
3. Adapter trimming
   Cutadapt
       │
       ▼
4. Reference alignment
   HISAT2
       │
       ▼
   SAM → sorted BAM
   SAMtools
       │
       ▼
5. Duplicate marking
   Sambamba
       │
       ▼
6. Mapping quality control
   SAMtools flagstat + stats
       │
       ▼
7. Gene-level quantification
   featureCounts / Rsubread
       │
       ▼
Raw gene-count matrix
```

---

# Analysis scope

This workflow ends with the generation of the **raw gene-count matrix**.

Downstream analyses such as:

* count normalisation;
* sample concordance assessment;
* differential-expression analysis;
* gene-set enrichment;

are performed separately.

No biological sample is excluded during generation of the raw count matrix. Any sample exclusion required for downstream statistical analyses should therefore be performed **after** this pipeline.

---

# Software

| Software                 | Version | Purpose                                  |
| ------------------------ | ------: | ---------------------------------------- |
| FastQC                   |  0.12.1 | Raw-read quality control                 |
| MultiQC                  |    1.14 | QC report aggregation                    |
| Cutadapt                 |     4.3 | Illumina adapter trimming                |
| HISAT2                   |   2.2.1 | RNA-seq read alignment                   |
| SAMtools                 |    1.19 | BAM conversion, sorting, indexing and QC |
| Sambamba                 |   1.0.1 | Duplicate marking                        |
| Rsubread / featureCounts |  2.16.1 | Gene-level read counting                 |

The examples below use **Singularity** containers. Apptainer can be used equivalently where available.

---

# Repository structure

```text
.
├── transcriptomics.csv.manifest
│
├── 0_refseq/
│   ├── reference.fa
│   └── annotation.gtf
│
├── 0_tools/
│
├── 0_samples/
│
├── 01_fastqc/
│
├── 02_multiqc/
│   └── results/
│
├── 03_cutadapt/
│
├── 04_hisat2/
│   ├── index/
│   └── bam/
│
├── 05_MarkDuplicates/
│   └── bam/
│
├── 05_quality/
│   ├── per_sample/
│   └── summary_metrics.csv
│
└── 11_featureCounts/
    ├── cnt.csv
    ├── featureCounts_assignment_summary.csv
    └── sessionInfo.txt
```

Commands below assume that they are run from the **repository root**.

---

# Input data

## Sample manifest

Samples are described in:

```text
transcriptomics.csv.manifest
```

using the following format:

```text
ID,R1,R2
Hc_T_R_ARA_F_1,ftp.sra.ebi.ac.uk/..._1.fastq.gz,ftp.sra.ebi.ac.uk/..._2.fastq.gz
Hc_T_R_ARA_F_2,ftp.sra.ebi.ac.uk/..._1.fastq.gz,ftp.sra.ebi.ac.uk/..._2.fastq.gz
```

The sample identifier provided in the `ID` column is retained throughout the workflow.

URLs without an explicit scheme can be interpreted as FTP addresses by prepending:

```text
ftp://
```

## Reference genome and annotation

Place the reference files in:

```text
0_refseq/
```

with one PRJEB506 WBPS18 genome:

```text
0_refseq/reference.fa
```

and the corresponding canonical gene annotation:

```text
0_refseq/annotation.gtf
```

The genome and GTF must originate from the same reference release.

---

# Container preparation

Create the tool directory:

```bash
ROOT="$(pwd -P)"

mkdir -p "${ROOT}/0_tools"
```

Download the version-pinned BioContainers used by the workflow:

```bash
ROOT="$(pwd -P)"

wget -c \
    -O "${ROOT}/0_tools/fastqc_0.12.1.sif" \
    "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"

wget -c \
    -O "${ROOT}/0_tools/multiqc_1.14.sif" \
    "https://depot.galaxyproject.org/singularity/multiqc%3A1.14--pyhdfd78af_0"

wget -c \
    -O "${ROOT}/0_tools/cutadapt_4.3.sif" \
    "https://depot.galaxyproject.org/singularity/cutadapt%3A4.3--py310h1425a21_0"

wget -c \
    -O "${ROOT}/0_tools/hisat2_2.2.1.sif" \
    "https://depot.galaxyproject.org/singularity/hisat2%3A2.2.1--h87f3376_4"

wget -c \
    -O "${ROOT}/0_tools/samtools_1.19.sif" \
    "https://depot.galaxyproject.org/singularity/samtools%3A1.19--h50ea8bc_0"

wget -c \
    -O "${ROOT}/0_tools/sambamba_1.0.1.sif" \
    "https://depot.galaxyproject.org/singularity/sambamba%3A1.0.1--h6f6fda4_0"

wget -c \
    -O "${ROOT}/0_tools/rsubread_2.16.1.sif" \
    "https://depot.galaxyproject.org/singularity/bioconductor-rsubread%3A2.16.1--r43ha9d7317_0"
```

These files only need to be downloaded once.

---

# 0 — Acquire sequencing reads

## Input

```text
transcriptomics.csv.manifest
```

## Output

```text
0_samples/<ID>/<ID>_R1.fastq.gz
0_samples/<ID>/<ID>_R2.fastq.gz
```

FASTQ files are downloaded using `wget -c`, allowing interrupted transfers to be resumed.

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

mkdir -p "${ROOT}/0_samples"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    r1="${r1//$'\r'/}"
    r2="${r2//$'\r'/}"

    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    case "${r1}" in
        *://*) ;;
        *) r1="ftp://${r1}" ;;
    esac

    case "${r2}" in
        *://*) ;;
        *) r2="ftp://${r2}" ;;
    esac

    mkdir -p "${ROOT}/0_samples/${id}"

    wget -c \
        -O "${ROOT}/0_samples/${id}/${id}_R1.fastq.gz" \
        "${r1}"

    wget -c \
        -O "${ROOT}/0_samples/${id}/${id}_R2.fastq.gz" \
        "${r2}"

    gzip -t "${ROOT}/0_samples/${id}/${id}_R1.fastq.gz"
    gzip -t "${ROOT}/0_samples/${id}/${id}_R2.fastq.gz"

done < "${MANIFEST}"
```

---

# 1 — Raw-read quality control

Raw paired-end reads are assessed with **FastQC 0.12.1** before trimming.

## Input

```text
0_samples/<ID>/<ID>_R1.fastq.gz
0_samples/<ID>/<ID>_R2.fastq.gz
```

## Output

```text
01_fastqc/<ID>/
```

containing the FastQC HTML and ZIP reports for both mates.

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

mkdir -p "${ROOT}/01_fastqc"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    mkdir -p "${ROOT}/01_fastqc/${id}"

    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/fastqc_0.12.1.sif" \
        fastqc \
        --threads 4 \
        --outdir "${ROOT}/01_fastqc/${id}" \
        "${ROOT}/0_samples/${id}/${id}_R1.fastq.gz" \
        "${ROOT}/0_samples/${id}/${id}_R2.fastq.gz"

done < "${MANIFEST}"
```

---

# 2 — Summarise raw-read QC

FastQC reports are aggregated using **MultiQC 1.14**.

## Input

```text
01_fastqc/
```

## Output

```text
02_multiqc/results/raw_reads_multiqc.html
```

```bash
ROOT="$(pwd -P)"

mkdir -p "${ROOT}/02_multiqc/results"

singularity exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/multiqc_1.14.sif" \
    multiqc \
    --filename raw_reads_multiqc.html \
    --outdir "${ROOT}/02_multiqc/results" \
    "${ROOT}/01_fastqc"
```

---

# 3 — Adapter trimming

Paired-end adapter trimming is performed with **Cutadapt 4.3**.

The Illumina adapter sequences are specified explicitly:

| Mate | Adapter                             |
| ---- | ----------------------------------- |
| R1   | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` |
| R2   | `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` |

Other Cutadapt trimming parameters are left at their default values.

## Input

```text
0_samples/<ID>/<ID>_R1.fastq.gz
0_samples/<ID>/<ID>_R2.fastq.gz
```

## Output

```text
03_cutadapt/<ID>/<ID>_R1.trim.fastq.gz
03_cutadapt/<ID>/<ID>_R2.trim.fastq.gz
```

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

mkdir -p "${ROOT}/03_cutadapt"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    mkdir -p "${ROOT}/03_cutadapt/${id}"

    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/cutadapt_4.3.sif" \
        cutadapt \
        -j 4 \
        -a "${ADAPTER_R1}" \
        -A "${ADAPTER_R2}" \
        -o "${ROOT}/03_cutadapt/${id}/${id}_R1.trim.fastq.gz" \
        -p "${ROOT}/03_cutadapt/${id}/${id}_R2.trim.fastq.gz" \
        "${ROOT}/0_samples/${id}/${id}_R1.fastq.gz" \
        "${ROOT}/0_samples/${id}/${id}_R2.fastq.gz"

    gzip -t "${ROOT}/03_cutadapt/${id}/${id}_R1.trim.fastq.gz"
    gzip -t "${ROOT}/03_cutadapt/${id}/${id}_R2.trim.fastq.gz"

done < "${MANIFEST}"
```

---

# 4 — Reference alignment

Trimmed reads are aligned to the *H. contortus* PRJEB506 WBPS18 genome using **HISAT2 2.2.1**.

Alignment output is subsequently converted to BAM, coordinate-sorted and indexed using **SAMtools 1.19**.

---

## 4.1 Build the HISAT2 index

## Input

```text
0_refseq/reference.fa
```

## Output

```text
04_hisat2/index/WBPS18/
```

```bash
ROOT="$(pwd -P)"
REFERENCE="$(find "${ROOT}/0_refseq" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fasta' \) -print -quit)"

mkdir -p "${ROOT}/04_hisat2/index/WBPS18"

singularity exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/hisat2_2.2.1.sif" \
    hisat2-build \
    -p 8 \
    "${REFERENCE}" \
    "${ROOT}/04_hisat2/index/WBPS18/hisat2idx"
```

---

## 4.2 Align paired-end reads

Read-group information is added during alignment.

For each sample:

```text
ID = sample ID
SM = sample ID
LB = library
PL = illumina
```

## Input

```text
03_cutadapt/<ID>/<ID>_R1.trim.fastq.gz
03_cutadapt/<ID>/<ID>_R2.trim.fastq.gz
```

## Output

```text
04_hisat2/bam/<ID>.sorted.bam
04_hisat2/bam/<ID>.sorted.bam.bai
```

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

mkdir -p "${ROOT}/04_hisat2/bam"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    SAM="${ROOT}/04_hisat2/bam/${id}.sam"
    BAM="${ROOT}/04_hisat2/bam/${id}.bam"
    SORTED="${ROOT}/04_hisat2/bam/${id}.sorted.bam"

    # Alignment
    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/hisat2_2.2.1.sif" \
        hisat2 \
        -p 8 \
        --rg-id "${id}" \
        --rg "SM:${id}" \
        --rg "LB:library" \
        --rg "PL:illumina" \
        -x "${ROOT}/04_hisat2/index/WBPS18/hisat2idx" \
        -1 "${ROOT}/03_cutadapt/${id}/${id}_R1.trim.fastq.gz" \
        -2 "${ROOT}/03_cutadapt/${id}/${id}_R2.trim.fastq.gz" \
        -S "${SAM}"

    # SAM → BAM
    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/samtools_1.19.sif" \
        samtools view \
        -@ 4 \
        -bS \
        -o "${BAM}" \
        "${SAM}"

    rm "${SAM}"

    # Coordinate sorting
    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/samtools_1.19.sif" \
        samtools sort \
        -@ 8 \
        -o "${SORTED}" \
        "${BAM}"

    rm "${BAM}"

    # BAM indexing
    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/samtools_1.19.sif" \
        samtools index \
        -@ 8 \
        "${SORTED}"

done < "${MANIFEST}"
```

HISAT2 alignment parameters other than threading, read-group information, index and paired-end inputs are left at their defaults in this workflow.

---

# 5 — Duplicate marking

Duplicates are marked using **Sambamba 1.0.1**.

Duplicates are **marked rather than removed**.

## Input

```text
04_hisat2/bam/<ID>.sorted.bam
```

## Output

```text
05_MarkDuplicates/bam/<ID>.mkdup.sorted.bam
05_MarkDuplicates/bam/<ID>.mkdup.sorted.bam.bai
```

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

mkdir -p "${ROOT}/05_MarkDuplicates/bam"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    BAM="${ROOT}/05_MarkDuplicates/bam/${id}.mkdup.sorted.bam"

    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/sambamba_1.0.1.sif" \
        sambamba markdup \
        -t=8 \
        "${ROOT}/04_hisat2/bam/${id}.sorted.bam" \
        "${BAM}"

    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/samtools_1.19.sif" \
        samtools index \
        -@ 8 \
        "${BAM}"

done < "${MANIFEST}"
```

---

# 6 — Mapping quality control

Mapping statistics are calculated from the duplicate-marked BAM files using **SAMtools 1.19**.

Both complete SAMtools reports and a compact cross-sample summary are retained.

## Input

```text
05_MarkDuplicates/bam/<ID>.mkdup.sorted.bam
```

## Per-sample outputs

```text
05_quality/per_sample/<ID>/<ID>_flagstat.txt
05_quality/per_sample/<ID>/<ID>_stats.txt
```

## Summary output

```text
05_quality/summary_metrics.csv
```

---

## 6.1 Generate SAMtools reports

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

mkdir -p "${ROOT}/05_quality/per_sample"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    mkdir -p "${ROOT}/05_quality/per_sample/${id}"

    BAM="${ROOT}/05_MarkDuplicates/bam/${id}.mkdup.sorted.bam"

    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/samtools_1.19.sif" \
        samtools flagstat \
        -@ 4 \
        "${BAM}" \
        > "${ROOT}/05_quality/per_sample/${id}/${id}_flagstat.txt"

    singularity exec --cleanenv \
        --bind "${ROOT}:${ROOT}" \
        "${ROOT}/0_tools/samtools_1.19.sif" \
        samtools stats \
        -@ 4 \
        "${BAM}" \
        > "${ROOT}/05_quality/per_sample/${id}/${id}_stats.txt"

done < "${MANIFEST}"
```

---

## 6.2 Generate a compact QC table

The summary contains:

| Column                | Description                                |
| --------------------- | ------------------------------------------ |
| `Sample`              | Sample identifier                          |
| `TotalReads`          | Total reads reported by `flagstat`         |
| `MappedReads`         | Mapped reads                               |
| `AlignmentRatePct`    | Percentage of mapped reads                 |
| `ProperlyPairedReads` | Properly paired reads                      |
| `ProperlyPairedPct`   | Properly paired percentage                 |
| `DuplicateReads`      | Reads marked as duplicates                 |
| `DuplicationRatePct`  | Duplicate reads / total reads × 100        |
| `BasesMappedCigar`    | Bases mapped according to the CIGAR string |

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

OUT="${ROOT}/05_quality/summary_metrics.csv"

printf 'Sample,TotalReads,MappedReads,AlignmentRatePct,ProperlyPairedReads,ProperlyPairedPct,DuplicateReads,DuplicationRatePct,BasesMappedCigar\n' > "${OUT}"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    FLAG="${ROOT}/05_quality/per_sample/${id}/${id}_flagstat.txt"
    STATS="${ROOT}/05_quality/per_sample/${id}/${id}_stats.txt"

    total="$(awk '/ in total / {print $1; exit}' "${FLAG}")"

    mapped_line="$(grep -m 1 ' mapped (' "${FLAG}")"
    proper_line="$(grep -m 1 ' properly paired (' "${FLAG}")"

    mapped="${mapped_line%% *}"
    proper="${proper_line%% *}"

    mapped_pct="$(sed -n 's/.*(\([0-9.][0-9.]*\)%.*/\1/p' <<< "${mapped_line}")"
    proper_pct="$(sed -n 's/.*(\([0-9.][0-9.]*\)%.*/\1/p' <<< "${proper_line}")"

    duplicates="$(awk '/ duplicates$/ {print $1; exit}' "${FLAG}")"

    duplicate_pct="$(
        awk \
            -v d="${duplicates:-0}" \
            -v t="${total:-0}" \
            'BEGIN {
                if (t > 0) printf "%.4f", 100*d/t;
                else print ""
            }'
    )"

    bases_cigar="$(
        awk -F '\t' \
            '$1=="SN" && $2 ~ /^bases mapped \(cigar\):/ {
                print $3;
                exit
            }' \
            "${STATS}"
    )"

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "${id}" \
        "${total}" \
        "${mapped}" \
        "${mapped_pct}" \
        "${proper}" \
        "${proper_pct}" \
        "${duplicates}" \
        "${duplicate_pct}" \
        "${bases_cigar}" \
        >> "${OUT}"

done < "${MANIFEST}"
```

---

# 7 — Gene-level quantification

Gene-level read counts are generated with **featureCounts** from **Rsubread 2.16.1**.

The analysis uses:

```text
paired-end counting
GTF annotation
gene-level count matrix
```

No strand-specific option is specified in this workflow.

Sample names in the resulting count matrix are taken directly from the manifest rather than inferred from BAM filenames.

## Input

```text
05_MarkDuplicates/bam/<ID>.mkdup.sorted.bam
0_refseq/annotation.gtf
```

## Main outputs

```text
11_featureCounts/cnt.csv
11_featureCounts/featureCounts_assignment_summary.csv
11_featureCounts/sessionInfo.txt
```

The count matrix has the form:

```text
GeneID,Length,<sample_1>,<sample_2>,...
HCON_00000010,...,...
HCON_00000020,...,...
...
```

---

## 7.1 Create the BAM manifest

```bash
ROOT="$(pwd -P)"
MANIFEST="${ROOT}/transcriptomics.csv.manifest"

mkdir -p "${ROOT}/11_featureCounts"

BAM_MANIFEST="${ROOT}/11_featureCounts/bam_manifest.tsv"

printf 'ID\tBAM\n' > "${BAM_MANIFEST}"

while IFS=, read -r id r1 r2; do

    id="${id//$'\r'/}"
    [[ "${id}" == "ID" || -z "${id}" ]] && continue

    printf '%s\t%s\n' \
        "${id}" \
        "${ROOT}/05_MarkDuplicates/bam/${id}.mkdup.sorted.bam" \
        >> "${BAM_MANIFEST}"

done < "${MANIFEST}"
```

---

## 7.2 Run featureCounts

Create:

```text
11_featureCounts/featureCounts.R
```

with:

```r
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L) {
    stop(
        "Usage: featureCounts.R ",
        "<annotation.gtf> ",
        "<bam_manifest.tsv> ",
        "<cnt.csv> ",
        "<assignment_summary.csv>"
    )
}

annotation   <- args[[1L]]
bam_manifest <- args[[2L]]
counts_out   <- args[[3L]]
summary_out  <- args[[4L]]

suppressPackageStartupMessages(
    library(Rsubread)
)

samples <- read.delim(
    bam_manifest,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

fc <- featureCounts(
    files = samples$BAM,
    annot.ext = annotation,
    isGTFAnnotationFile = TRUE,
    isPairedEnd = TRUE,
    nthreads = 4
)

colnames(fc$counts) <- samples$ID

count_table <- data.frame(
    fc$annotation[, c("GeneID", "Length"), drop = FALSE],
    fc$counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

write.table(
    count_table,
    file = counts_out,
    quote = FALSE,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)

assignment <- as.data.frame(
    fc$stat,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

colnames(assignment)[-1L] <- samples$ID

write.table(
    assignment,
    file = summary_out,
    quote = FALSE,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE
)

writeLines(
    c(
        paste0("Rsubread=", packageVersion("Rsubread")),
        capture.output(sessionInfo())
    ),
    con = file.path(
        dirname(counts_out),
        "sessionInfo.txt"
    )
)
```

Run it with:

```bash
ROOT="$(pwd -P)"

ANNOTATION="$(find "${ROOT}/0_refseq" -maxdepth 1 -type f -name '*.gtf' -print -quit)"

singularity exec --cleanenv \
    --bind "${ROOT}:${ROOT}" \
    "${ROOT}/0_tools/rsubread_2.16.1.sif" \
    Rscript \
    "${ROOT}/11_featureCounts/featureCounts.R" \
    "${ANNOTATION}" \
    "${ROOT}/11_featureCounts/bam_manifest.tsv" \
    "${ROOT}/11_featureCounts/cnt.csv" \
    "${ROOT}/11_featureCounts/featureCounts_assignment_summary.csv"
```

---

# Important methodological points

## Adapter trimming

Only the two Illumina adapter sequences are explicitly provided to Cutadapt.

```text
Cutadapt 4.3
   ├── R1 Illumina adapter
   ├── R2 Illumina adapter
   └── otherwise default trimming behaviour
```

No additional quality or minimum-length filtering parameter is introduced at this stage.

## Alignment

HISAT2 uses the PRJEB506 WBPS18 reference genome.

```text
trimmed paired-end FASTQ
          │
          ▼
       HISAT2
          │
          ▼
         SAM
          │
          ▼
     SAMtools view
          │
          ▼
     SAMtools sort
          │
          ▼
   coordinate-sorted BAM
```

## Duplicate handling

Sambamba is used to **mark duplicates**.

The resulting duplicate-marked BAM files are then used both for mapping QC and featureCounts.

```text
sorted BAM
    │
    ▼
Sambamba markdup
    │
    ├── SAMtools QC
    │
    └── featureCounts
```

Duplicates are not explicitly removed before gene counting in this workflow.

## Counting

featureCounts is run in paired-end mode:

```r
isPairedEnd = TRUE
```

with the WBPS18 GTF annotation.

Other featureCounts parameters are not explicitly modified and therefore follow Rsubread defaults.

---

# Main outputs

| Analysis                            | Main output                                             |
| ----------------------------------- | ------------------------------------------------------- |
| Raw-read QC                         | `01_fastqc/<ID>/`                                       |
| Aggregated raw-read QC              | `02_multiqc/results/raw_reads_multiqc.html`             |
| Trimmed reads                       | `03_cutadapt/<ID>/<ID>_R1.trim.fastq.gz`                |
| Trimmed reads                       | `03_cutadapt/<ID>/<ID>_R2.trim.fastq.gz`                |
| Aligned reads                       | `04_hisat2/bam/<ID>.sorted.bam`                         |
| Duplicate-marked reads              | `05_MarkDuplicates/bam/<ID>.mkdup.sorted.bam`           |
| Mapping QC                          | `05_quality/per_sample/<ID>/`                           |
| Mapping summary                     | `05_quality/summary_metrics.csv`                        |
| Raw gene counts                     | `11_featureCounts/cnt.csv`                              |
| featureCounts assignment statistics | `11_featureCounts/featureCounts_assignment_summary.csv` |
| R environment information           | `11_featureCounts/sessionInfo.txt`                      |

---

# Downstream analysis

The principal output for downstream transcriptomic analysis is:

```text
11_featureCounts/cnt.csv
```

This table contains **raw, unnormalised gene counts**.

Any subsequent operations such as:

```text
sample QC
    ↓
sample exclusion where justified
    ↓
normalisation
    ↓
differential-expression modelling
    ↓
functional / gene-set enrichment
```

belong to the downstream statistical analysis and are deliberately kept separate from this preprocessing and quantification workflow.

---

# HPC execution

The commands documented here describe the **scientific analysis** rather than the complete production scheduler logic.

In the original HPC implementation:

* individual samples are processed as separate LSF jobs;
* independent samples can be run in parallel;
* temporary directories are used during execution;
* expected outputs are checked before a stage is accepted;
* failed LSF jobs prevent progression to subsequent stages;
* existing completed stages can be skipped.

Those mechanisms improve robustness on a shared HPC system but are intentionally omitted here to keep the analysis workflow easy to read and adapt.
