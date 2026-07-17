# Documentation of the transcriptomic pipeline and differential expression analysis

> **Purpose of the document.** This page converts a raw laboratory notebook into scientific documentation intended to accompany the code repository associated with a publication. It describes the **steps, commands, parameters, tools, versions, files, conventions and activation states recorded in the supplied snapshot**.  
>
> The purpose is not to present the repository as an immediately re-runnable workflow, but to allow readers to understand precisely **what was coded and how the analyses were conducted or prepared**.

## General metadata

| Field | Recorded information |
|---|---|
| Organism | *Haemonchus contortus* |
| Data | Paired-end Illumina RNA-seq |
| Main analysis | genome alignment, gene-level counting and differential expression |
| Reference genome | `haemonchus_contortus.PRJEB506.WBPS18.genomic.fa` |
| Main annotation | `haemonchus_contortus.PRJEB506.WBPS18.canonical_geneset.gtf` |
| Alternative transcriptome | `haemonchus_contortus.PRJEB506.WBPS18.mRNA_transcripts.fa` |
| Main infrastructure | HPC cluster with a Slurm environment and software modules |
| Alternative infrastructure | Singularity images for FastQC, MultiQC, Cutadapt, HISAT2 and samtools |
| R version invoked | R 4.3.2 |
| Author named in the historical README | Robin Lioutaud |
| Unit named | UMR INRAE/ENVT 1436 InTheRes, Toulouse, France |
| Licence stated | CeCILL |
| Associated publication | To be provided |
| DOI / read accession | To be provided |
| Exact analysis date | Not recorded |
| Commit corresponding to the published results | To be provided |

## Scope of this compilation

The notebook contains an RNA-seq pipeline extending from quality control to differential expression visualisation. It includes:

- a main route using **FastQC → MultiQC → Cutadapt → HISAT2/samtools → Sambamba → featureCounts → DESeq2 → volcano plots**;
- an alternative pseudo-alignment branch using **Kallisto → Sleuth**;
- an alternative **edgeR** analysis;
- several exploratory normalisations or transformations: **CPM, MRM/DESeq2, TMM, TPM and RPKM**;
- sample exploration using **PCA** and a distance heatmap;
- archived tests with **HTSeq** and **MMQuant**.

No functional enrichment, GSEA, variant-calling or annotation-database integration script is present in this compilation, although the historical README mentions a ‘manual GSEA’.

### Interpretation of conditional blocks

The scripts were developed iteratively:

- `if (TRUE)` or `if true`: active block in the snapshot;
- `if (FALSE)` or `if false`: disabled block in the snapshot;
- a disabled block may correspond to a step that had been run previously when its outputs are subsequently read;
- the state of the code alone therefore does not demonstrate that every step was executed for each result in the publication.

## Overview

```text
Raw paired-end FASTQ files
        │
        ├── FastQC ── MultiQC
        │
        └── Cutadapt
                │
                └── HISAT2 + samtools
                        │
                        └── sorted and indexed BAM files
                                │
                                └── Sambamba markdup
                                        │
                                        └── featureCounts
                                                │
                                                ├── CPM / MRM / TMM / TPM
                                                ├── PCA + distance heatmap
                                                └── DESeq2
                                                        │
                                                        └── Plotly volcano plots

Alternative branch:
FASTQ → Kallisto → transcript abundances → Sleuth
```

## Documented analytical directory structure

```text
.
├── 01_fastqc/
├── 02_multiqc/
├── 03_cutadapt/
├── 04_hisat2/
│   └── autre/kallisto/
├── 05_MarkDuplicates/
├── 11_featurecounts/
│   └── autres/{htseq,mmquant}/
├── 13_normalisation_cpm/
├── 13_normalisation_mrm/
├── 13_normalisation_tmm/
├── 13_normalisation_tpm/
├── 14_normalisation/
├── 15_exploration/
├── 16_deseq2/
│   └── autres/{edgeR,sleuth}/
├── 17_DGEplot/
├── main.sh
└── README.md
```

The `.git/` directory was also present in the raw compilation. Its internal files and the example hooks supplied by Git are not analytical steps and are not reproduced in the appendix.

## Step status in the snapshot

| Directory | Function | Visible status |
|---|---|---|
| `01_fastqc` | initial quality control | active |
| `02_multiqc` | report aggregation | active |
| `03_cutadapt` | adapter removal | active |
| `04_hisat2` | index construction, alignment, sorting, BAM indexing and metrics | active |
| `04_hisat2/autre/kallisto` | alternative quantification | indexing and quantification disabled; TSV concatenation active |
| `05_MarkDuplicates` | duplicate marking | active |
| `05_MarkDuplicates/quality.sh` | post-duplicate metrics | `flagstat/stats` generation disabled; CSV summary active |
| `11_featurecounts` | gene-level counting | featureCounts calculation disabled; renaming of an existing `cnt.csv` active |
| `13_normalisation_*` | CPM, MRM, TMM and TPM | active when each script is run |
| `14_normalisation` | RPKM draft and normalisation selection | RPKM calculation active in the draft; launcher configured for `tpm.r`, which is absent from this compilation |
| `15_exploration` | PCA and distance heatmap | function active; input-file selection requires verification |
| `16_deseq2` | main differential expression analysis | active |
| `16_deseq2/autres/edgeR` | alternative multifactorial analysis | stand-alone script |
| `16_deseq2/autres/sleuth` | exploration of Kallisto results | stand-alone script |
| `17_DGEplot` | volcano plots | active; DGE heatmap disabled |
| `main.sh` | general orchestration | prepared but incompatible with the current directory names |

# Data organisation

The historical README specifies two sibling directories:

```text
<root>/
├── RAW/   # FASTQ files, adapters, genome and annotations
└── RES/   # pipeline repository and results by step
```

## Explicitly referenced inputs

| Category | Names or patterns used |
|---|---|
| Raw FASTQ files | `../../RAW/*.fastq*` |
| Absolute FASTQ path | `/work/user/rlioutaud/RAW/hcon/rnaseq` |
| Read pairs | `<sample>_1.fastq` and `<sample>_2.fastq` |
| Adapters | `../../RAW/adapter1*.txt`, `../../RAW/adapter2*.txt` |
| Genome | `*.genomic.fa*` |
| FASTA index | `haemonchus_contortus.PRJEB506.WBPS18.genomic.fa.fai` |
| GTF annotation | `haemonchus_contortus.PRJEB506.WBPS18.canonical_geneset.gtf` |
| GFF3 annotation | `haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3` |
| Kallisto transcriptome | `haemonchus_contortus.PRJEB506.WBPS18.mRNA_transcripts.fa` |
| Lane-to-sample mapping | `lane_sample_IDs.list` |
| Kallisto/Sleuth metadata | `rnaseq_metadata.txt` |
| Annotation-to-function mapping | `AnnotationToGeneFunction.csv` |
| Main counts | `11_featurecounts/cnt.csv` |

## Sample-naming convention

The code uses dot-delimited names, for example:

```text
T.R.ARA.F.1
```

| Position or code | Observed values | Interpretation used |
|---|---|---|
| `R` / `S` | resistance status | Resistant / Susceptible |
| `ARA`, `BET`, `BUN`, `MOU`, `CHI`, `LUC` | population or strain | biological groups |
| `F` / `M` | sex | Female / Male |
| terminal integer | replicate | replicate number |
| first element, for example `T` | undefined | meaning not specified |

The README gives `R.LOC.M.1.fq.gz` as an example, whereas the preprocessing scripts mainly use the `_1.fastq` and `_2.fastq` suffixes.

# Software, environments and versions

| Software / environment | Recorded version | Invocation mode | Note |
|---|---:|---|---|
| Ubuntu | 20.04 | base image for the Singularity definitions | FastQC, MultiQC, Cutadapt and HISAT2 containers |
| FastQC | 0.12.1 | HPC module | container: 0.11.9 |
| MultiQC | 1.14 | HPC module | unpinned `pip3 install multiqc` installation in the container |
| Cutadapt | 4.3 | HPC module | unpinned `pip3 install cutadapt` installation in the container |
| HISAT2 | 2.2.1 | HPC module | container: 2.2.0 |
| samtools | 1.19 | HPC module | used for conversion, sorting, indexing and metrics |
| Python | 3.11.1 | HPC module | loaded during the HISAT2 step |
| Kallisto | 0.50.1 | HPC module | alternative branch |
| Sambamba | 1.0.1 | HPC module | duplicate marking |
| MMQuant | 1.0.9 | HPC module | alternative counting method |
| R | 4.3.2 | absolute executable path | `/save/user/...` and `/work/user/...` paths |
| Singularity | not recorded | `.sif` images | engine version absent from the notebook |
| Rsubread / featureCounts | not recorded | R package | package version absent |
| DESeq2 | not recorded | R package | main analysis and MRM normalisation |
| edgeR | not recorded | R package | TMM and alternative analysis |
| Sleuth | not recorded | R package | Kallisto branch |
| ggplot2 | not recorded | R package | plotting |
| Plotly | not recorded | R package | interactive PCA and volcano plots |
| pheatmap | not recorded | R package | heatmaps |
| tidyverse | not recorded | R package | data manipulation |
| RColorBrewer | not recorded | R package | heatmap palettes |
| patchwork | not recorded | R package | assembly of Sleuth PCA plots |

No `sessionInfo()`, `renv.lock` file, Conda environment or frozen Python dependency specification is present.

# Main pipeline

## 1. Initial quality control — FastQC

**Directory:** `01_fastqc/`  
**Inputs:** paired-end FASTQ files in `../../RAW/`  
**Outputs:** FastQC reports in the current directory.

### HPC

```bash
module load bioinfo/FastQC/0.12.1

fastqc -t 4 -o ./ ../../RAW/${file}_1.fastq* &
fastqc -t 4 -o ./ ../../RAW/${file}_2.fastq*
```

Identifiers are generated by removing everything after the first underscore:

```bash
FILENAME=$(echo "$FILENAME" | sed 's/_.*//')
```

Samples are divided into groups of six. For each sample, both mates are processed within a background subprocess. The theoretical maximum is therefore twelve FastQC processes using four threads each, i.e. 48 threads.

### Singularity

```bash
singularity exec ../fastqc.sif fastqc -t 4 -o ./ ../../RAW/${file}_1.fastq*
singularity exec ../fastqc.sif fastqc -t 4 -o ./ ../../RAW/${file}_2.fastq*
```

The container definition installs FastQC 0.11.9 on Ubuntu 20.04 with Java 11.

## 2. Quality-control aggregation — MultiQC

**Directory:** `02_multiqc/`

```bash
module load bioinfo/MultiQC/1.14
multiqc --outdir ./ ../*fastqc
```

Alternative:

```bash
singularity exec ../multiqc.sif multiqc --outdir ./ ../*fastqc
```

The version installed by `pip3 install multiqc` in the Singularity definition is not pinned.

## 3. Adapter removal — Cutadapt

**Directory:** `03_cutadapt/`  
**Inputs:** paired-end FASTQ and adapter files  
**Outputs:** `<sample>_1.trim.fastq`, `<sample>_2.trim.fastq`.

```bash
module load bioinfo/Cutadapt/4.3

adapter1=$(cat ../../RAW/adapter1*.txt | tr -d '\n' | tr -d '\r')
adapter2=$(cat ../../RAW/adapter2*.txt | tr -d '\n' | tr -d '\r')

cutadapt \
  -a $adapter1 \
  -A $adapter2 \
  -o ${samplename}_1.trim.fastq \
  -p ${samplename}_2.trim.fastq \
  $inputpath/${fastq}_1.fastq \
  $inputpath/${fastq}_2.fastq
```

No explicit quality threshold, minimum length, error rate or maximum number of ambiguous bases is specified. On the HPC system, all pairs are launched in the background without an explicit concurrency limit, despite `#SBATCH --ntasks=1`.

## 4. Alignment — HISAT2 and samtools

**Directory:** `04_hisat2/`  
**Inputs:** genome and trimmed FASTQ files  
**Outputs:** temporary SAM files, sorted BAM files, BAM indices and metrics.

### Path detection

```bash
PathGenome=$(find "../../RAW/" -name "*.genomic.fa*" -print -quit)
Path=$(find ../../RES -type d -name "*cutadapt*" -print)
```

If several directories match `*cutadapt*`, the `Path` variable may contain multiple lines.

### Index construction

```bash
IndexName="hisat2idx"
hisat2-build -p 8 $PathGenome $IndexName
chmod 777 hisat2idx*
```

### Paired-end alignment and read groups

```bash
hisat2 \
  -p 8 \
  --rg ID:"$fastq" \
  --rg SM:"$fastq" \
  --rg LB:"library" \
  --rg PL:"illumina" \
  -x hisat2idx \
  -1 ${Path}/${fastq}_1.trim.fastq \
  -2 ${Path}/${fastq}_2.trim.fastq \
  -S "./${fastq}.sam"
```

No library-strandedness or known-transcriptome option (`--rna-strandness`, splice-site file or exon file) is recorded.

### Conversion, sorting, indexing and metrics

```bash
samtools view -bS "${fastq}.sam" |
  samtools sort -@ 4 -o "${fastq}.sorted.bam"

rm "${fastq}.sam"
samtools index -@ 4 "$bam"
samtools flagstat "$BAM_FILE" > "${BASENAME}_flagstat.txt"
samtools stats "$BAM_FILE" > "${BASENAME}_stats.txt"
```

Four samples are launched simultaneously and each HISAT2 process receives eight threads. The theoretical requirement therefore reaches 32 threads, whereas the Slurm comment specifies `--ntasks=8`.

## 5. Duplicate marking — Sambamba

**Directory:** `05_MarkDuplicates/`  
**Inputs:** `*.sorted.bam`  
**Outputs:** `*.mkdup.sorted.bam`.

```bash
module load bioinfo/Sambamba/1.0.1

sambamba markdup -t=8 \
  "${input}/${bam}.sorted.bam" \
  "./${bam}.mkdup.sorted.bam"
```

Three BAM files are processed in parallel, corresponding to up to 24 threads. Duplicates are **marked**, not removed. No indexing of the post-marking BAM files is explicitly recorded.

### Quality summary

The following blocks are disabled:

```bash
samtools flagstat "$BAM_FILE" > "${BASENAME}_flagstat.txt"
samtools stats "$BAM_FILE" > "${BASENAME}_stats.txt"
```

The active block nevertheless reads these files and writes `summary_metrics.csv`:

```text
File,Total Reads,Mapped Reads,Alignment Rate,
Uniquely Mapped Reads,Duplication Rate,Coverage Depth
```

Main calculations:

```bash
alignment_rate = mapped_reads / total_reads * 100
uniquely_mapped_reads = (mapped_reads - secondary) / total_reads * 100
duplication_rate = duplicates / total_reads * 100
coverage_depth = bases_mapped_cigar / genome_size
```

The term `Uniquely Mapped Reads` is an internal approximation based on subtracting secondary alignments rather than the ‘uniquely mapped’ metric produced by an aligner. Variables named `rrna_rate` and `exonic_rate` are calculated from `QC-passed reads` and `QC-failed reads` lines, but are not exported and do not directly measure rRNA content or exonic mapping.

## 6. Gene-level counting — featureCounts

**Directory:** `11_featurecounts/`  
**Expected input:** post-marking BAM files  
**Output:** `cnt.csv`.

```r
bams <- list.files(
  "../../RES",
  pattern = ".*MarkDuplicates.*\\.bam$",
  full.names = TRUE,
  recursive = TRUE
)

counts <- featureCounts(
  files = bams,
  annot.ext = "../../RAW/haemonchus_contortus.PRJEB506.WBPS18.canonical_geneset.gtf",
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  nthreads = 4
)
```

Export:

```r
write.table(
  x = data.frame(
    counts$annotation[, c("GeneID", "Length")],
    counts$counts,
    stringsAsFactors = FALSE
  ),
  file = "cnt.csv",
  quote = FALSE,
  sep = ",",
  row.names = FALSE
)
```

The calculation is disabled (`if (FALSE)`). The active block re-reads an existing `cnt.csv` and shortens the sample names by removing their common prefix and suffix.

The `.*MarkDuplicates.*\\.bam$` pattern requires verification: if `MarkDuplicates` occurs only in the directory name and not in the BAM filename, the search may select no files.

No explicit library-strandedness, feature-type or GTF-attribute parameter is supplied; the default `featureCounts` values apply to omitted options.

### Archived alternatives

#### HTSeq

Fully commented example command:

```bash
htseq-count -f bam -r pos -s no -i gene_id output.bam reference.gtf > counts.txt
```

#### MMQuant 1.0.9

```bash
module load bioinfo/mmquant/1.0.9
mmquant --bam "$bam" --annotation "$annotation" --output "${filename}_counts.csv"
```

The paths in this branch (`../../2xAlignement/hisat2/`, `../featureCounts/annotations.gtf`) belong to an older directory organisation.

## 7. Alternative quantification — Kallisto

**Directory:** `04_hisat2/autre/kallisto/`

```bash
module load bioinfo/kallisto/0.50.1
```

### Index

```bash
kallisto index --index TRANSCRIPTS.ixd $PathTranscriptome
```

### Quantification

```bash
kallisto quant \
  --index ./TRANSCRIPTS.ixd \
  --output-dir kallisto_${SAMPLE_ID}_out \
  --bootstrap-samples 100 \
  --threads 7 \
  ${PathRnaseq}/${LANE}_1*.fastq \
  ${PathRnaseq}/${LANE}_2*.fastq \
  >> kallisto_${SAMPLE_ID}.log 2>&1
```

Index construction and quantification are disabled in the snapshot. The active phase extracts the fourth column of the `abundance.tsv` files, normally corresponding to `est_counts`, and then pastes the columns into `counts.tsv`:

```bash
cut -f4 -d$'\t' "$tsv" | tail -n +2 > "${filename}.part.tsv"
paste *.part.tsv > counts.tsv
```

Transcript identifiers are not retained in this concatenation. An identical row order across all files is therefore implicitly assumed.

## 8. Normalisations and exploratory tables

All methods start from a `cnt.csv` file containing `GeneID`, `Length` and the counts.

| Directory | Method | Output | Formula or function |
|---|---|---|---|
| `13_normalisation_cpm` | CPM | `cnt.cpm.csv` | `count / total_mapped × 10^6` |
| `13_normalisation_mrm` | DESeq2 size factors | `cnt.mrm.csv` | `estimateSizeFactors()` followed by `counts(normalized=TRUE)` |
| `13_normalisation_tmm` | ‘de-CPMised’ edgeR TMM | `cnt.tmm.csv` | TMM-CPM followed by multiplication by `TotalMapped / 10^6` |
| `13_normalisation_tpm` | manual TPM | `cnt.tpm.csv` | division by length followed by rescaling the sum to `10^6` |
| `14_normalisation` | RPKM draft | `rpkm.csv` or `tmm_rpkm.csv` | CPM × `1000 / Length` |

The MRM normalisation creates a DESeq2 object with a `~1` design solely to estimate size factors. The main differential expression analysis starts again from the **raw counts** and estimates its own factors.

The `14_normalisation/sh.sh` script calls `tpm.r`, although this file is absent from the compilation. Calls to `tmm.r`, `cpm.r` and `mrm.r` are commented out.

## 9. Sample exploration — PCA and heatmap

**Directory:** `15_exploration/`  
**Intended outputs:** `plots.html`, `ClusteredHeatmap.svg`.

Steps:

1. read a normalised table;
2. remove the length column;
3. transpose to samples × genes;
4. add `epsilon = 1`;
5. apply an optional transformation;
6. remove invariant genes;
7. run `prcomp(..., scale. = TRUE)`;
8. generate Plotly PCA plots coloured according to each component of the sample name;
9. calculate `1 - cor` distances and generate a hierarchical heatmap.

Coded transformations:

```r
if (tf == "rlog") {
  df <- log(df)
} else if (tf == "vst") {
  variances <- apply(df, 2, var)
  vst_transform <- sqrt(variances)
  df <- log(df + 1) / vst_transform
}
```

These transformations are **not** the DESeq2 `rlog()` and `vst()` functions.

The percentage of variance is recalculated as follows:

```r
varexp <- summary(pca)$importance[2, ]
varexp <- varexp^2
varexp <- varexp / sum(varexp) * 100
```

The variance proportion supplied by `prcomp` is therefore squared a second time; the displayed values do not correspond to conventional variance proportions.

The final call constructs a literal path:

```r
file.path(
  list.files("../", pattern = ".*normalisation_mrm.*", full.names = TRUE),
  ".*mrm.csv"
)
```

It does not actually search for files using a regular expression within the directory. This point should be checked against the historical execution.

## 10. Main differential expression analysis — DESeq2

**Directory:** `16_deseq2/`  
**Input:** `../11_featurecounts/cnt.csv`  
**Main output:** `cnt.deseq2.dge.csv`.

### Preparation

```r
counts <- read.csv(file, header = TRUE, row.names = 1)
counts <- counts[, -1]
counts <- as.matrix(counts)
counts <- matrix(
  as.integer(counts),
  nrow = nrow(counts),
  ncol = ncol(counts),
  dimnames = list(rownames(counts), colnames(counts))
)
```

### Selected contrast

The active code retains females:

```r
counts <- SeleCols(counts, "\\.F\\.")
```

Reference:

```r
reference <- SeleCols(reference, "\\.CHI\\.|\\.LUC\\.")
```

Comparison:

```r
comparison <- SeleCols(comparison, "\\.BET\\.")
```

The contrast is therefore:

```text
BET females versus CHI/LUC females
```

provided that the column names follow the expected conventions exactly.

### Model and pre-filtering

```r
conditions <- factor(
  c(rep("ref", ncol(reference)), rep("comp", ncol(comparison))),
  levels = c("ref", "comp")
)

dds <- DESeqDataSetFromMatrix(
  countData = CountsRefComp,
  colData = data.frame(condition = conditions),
  design = ~ condition
)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]
```

### Normalisation and testing

```r
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]

write.csv(
  res,
  file = paste0(filename, ".deseq2.dge.csv"),
  row.names = TRUE
)
```

With the levels ordered as `ref` followed by `comp`, the default coefficient normally represents `comp` relative to `ref`. A positive `log2FoldChange` therefore corresponds to higher expression in BET females than in CHI/LUC females.

The script temporarily exports normalised counts to `11_featurecounts/cnt.mrm.csv`, then removes all files in the source directory whose names match `.mrm.csv`.

### Dispersion plots

```r
pdf("cnt_dispersion.pdf")
plotDispEsts(dds)
dev.off()
```

A second block attempts to produce a custom SVG from `dispersionFunction(dds)`. This object is a dispersion function, and the resulting graph should not be equated with the standard `plotDispEsts()` output without verification.

The MA plot is present but disabled.

## 11. Alternative analysis — edgeR

**Directory:** `16_deseq2/autres/edgeR/`

The script reads an older table:

```r
/save/user/rlioutaud/pipelines/rnaseq/3vCounting/featureCounts/counts.txt
```

Positions 2 to 4 of the sample names are converted into factors, after which the following model is constructed:

```r
dge <- DGEList(counts = counts, group = coldata$RS)
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + RS + FM + RS:FM, data = coldata)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
results <- topTags(qlf)
```

Points to note:

- the initial comment mentions DESeq2, whereas the package used is edgeR;
- no explicit `estimateDisp()` step is present before `glmQLFit()`;
- `coef=2` refers to the second column of the design matrix without explicitly naming the contrast;
- only the results printed in `log.log` are recorded; no CSV export is performed.

This branch should be regarded as exploratory.

## 12. Kallisto/Sleuth branch

**Directory:** `16_deseq2/autres/sleuth/`

```r
metadata <- read.table("rnaseq_metadata.txt", header = TRUE)
so <- sleuth_prep(
  metadata,
  extra_bootstrap_summary = TRUE,
  num_cores = 2
)
```

The script:

- extracts an `est_counts` matrix;
- performs PCA;
- joins the scores to the covariates;
- generates three PCA plots coloured by `population`, `sex` and `drug`;
- calculates a Jensen–Shannon divergence matrix from TPM values;
- produces an annotated heatmap.

Outputs:

```text
figure_rnaseq_pca_pop_sex_drug.png
figure_rnaseq_heatmap_pop_sex_drug.pdf
figure_rnaseq_heatmap_pop_sex_drug.png
```

No Sleuth differential model (`sleuth_fit`, `sleuth_wt` or `sleuth_lrt`) is recorded. This branch therefore primarily performs sample exploration.

## 13. Volcano plots — Plotly

**Directory:** `17_DGEplot/`  
**Inputs:** `../16_deseq2/*.dge.csv`  
**Outputs:** `*_volcano_plot.html`.

```r
names(dge)[names(dge) == "log2FoldChange"] <- "log2FC"
names(dge)[names(dge) == "padj"] <- "neglog10padj"
dge$neglog10padj <- -log10(dge$neglog10padj)
```

Thresholds:

| Class | Criterion |
|---|---|
| non-significant or small effect | `abs(log2FC) < 1` or `-log10(padj) < 1.3` |
| downregulated | `log2FC <= -1` and `-log10(padj) >= 1.3` |
| upregulated | `log2FC >= 1` and `-log10(padj) >= 1.3` |

`1.3` corresponds approximately to `padj = 0.05`.

```r
htmlwidgets::saveWidget(
  volcan,
  paste0(filename, "_volcano_plot.html")
)
```

The DGE heatmap is disabled. It also depends on `dds` and `res` objects that are not created in this script and on the external files `Mcounts_normalized.csv` and `Fcounts_normalized.csv`.

# Output inventory

| Step | Main outputs |
|---|---|
| FastQC | HTML reports and ZIP archives |
| MultiQC | MultiQC report |
| Cutadapt | `*_1.trim.fastq`, `*_2.trim.fastq` |
| HISAT2/samtools | `*.sorted.bam`, `.bai`, `*_flagstat.txt`, `*_stats.txt` |
| Sambamba | `*.mkdup.sorted.bam` |
| Post-duplicate quality | `summary_metrics.csv` |
| featureCounts | `cnt.csv` |
| CPM | `cnt.cpm.csv` |
| MRM | `cnt.mrm.csv` |
| TMM | `cnt.tmm.csv` |
| TPM | `cnt.tpm.csv` |
| Exploratory RPKM | `rpkm.csv`, `tmm_rpkm.csv` |
| PCA/heatmap | `plots.html`, `ClusteredHeatmap.svg` |
| DESeq2 | `cnt.deseq2.dge.csv`, `cnt_dispersion.pdf`, `cnt_dispersion.svg`, `log.log` |
| Volcano plot | `cnt.deseq2.dge_volcano_plot.html` |
| Kallisto | `abundance.tsv`, `abundance.h5`, `run_info.json`, `counts.tsv` |
| Sleuth | PCA PNG, heatmap PDF and PNG |

# Points requiring caution

1. **Different versions on the HPC system and in containers.** FastQC 0.12.1/0.11.9 and HISAT2 2.2.1/2.2.0.
2. **Unpinned dependencies.** MultiQC and Cutadapt are installed with `pip` without version constraints; R package versions are not recorded.
3. **Absolute paths.** Several scripts are tied to the `/work/user/rlioutaud` and `/save/user/rlioutaud` accounts.
4. **Misplaced Slurm directives.** The `#SBATCH` lines occur after a shell instruction; Slurm normally interprets them only when they precede every executable instruction.
5. **Parallelisation inconsistent with requested resources.** HISAT2 may use 32 threads for a commented request of eight tasks; Cutadapt launches all samples simultaneously.
6. **Very permissive file permissions.** `chmod 777` is applied to numerous files.
7. **Identifiers truncated at the first underscore.** `sed 's/_.*//'` may merge or truncate complex identifiers.
8. **featureCounts is not recalculated in the snapshot.** The active script assumes that `cnt.csv` is already present.
9. **BAM selection requires verification.** The featureCounts pattern includes `MarkDuplicates` in the filename pattern.
10. **Library strandedness is absent.** No library orientation is recorded for HISAT2 or featureCounts.
11. **Post-duplicate QC is inconsistent in the current state.** The active summary depends on metrics whose generation is disabled.
12. **Kallisto concatenation lacks identifiers.** `counts.tsv` relies solely on row order.
13. **Non-standard transformations labelled rlog/vst.** They do not correspond to the DESeq2 transformations.
14. **Non-standard PCA variance percentage.** The proportions are squared twice.
15. **Hard-coded DESeq2 contrast.** Only BET females versus CHI/LUC females are tested.
16. **Incomplete edgeR output.** No explicit dispersion estimation or final CSV is present.
17. **DGE heatmap is not self-contained.** The disabled block depends on external objects and files.
18. **General orchestration is not functional as written.** `main.sh` does not target the `01_fastqc`, `02_multiqc`, etc. directories.
19. **GSEA is announced but absent.** No GSEA script is present in this compilation.
20. **The exact date and provenance of the results are not recorded.** The snapshot must be linked manually to the commit and publication figures.

# Overall execution and Git repository

The `main.sh` script searches for directories using:

```bash
if [[ $folder =~ ^[0-9]*v ]]; then
```

The actual directories are named `01_fastqc`, `02_multiqc`, etc. They do not contain the letter `v` after their numeric prefix. The script therefore does not select them in its current state.

The historical README explicitly states that the pipeline was not finalised and that the steps had to be run manually. It also mentions a risk of Windows carriage returns:

```bash
sed -i 's/\r//' LeScriptProblematique.sh
```

# Information to complete before publication

| Item | Status |
|---|---|
| Article DOI | to be provided |
| SRA/ENA read accession | to be provided |
| Sample table | to be added |
| Execution date or period | to be added if available |
| Git commit associated with the results | to be fixed |
| Exact contrast presented in the article | to be confirmed |
| Files used to generate each figure | to be listed |
| FastQC/MultiQC reports | to be archived or linked |
| R package versions | to be supplied only if they can be recovered historically |
| Reference checksums | to be added if available |

The historical scripts should preferably be retained unchanged in a `scripts/` subdirectory, with this document used as an interpreted description. Inconsistencies should not be silently corrected in the documentation of an analysis that has already been conducted.

# Appendix — analytical files reproduced without modification

This appendix reproduces the analytical files, launchers, container definitions, `main.sh` and the historical README. Internal `.git/` files are excluded.


## `01_fastqc/fastqc.def`

Source brute : lignes 5–42 de la compilation.

```text
Bootstrap: docker
From: ubuntu:20.04

%post
    export DEBIAN_FRONTEND=noninteractive
    echo "tzdata tzdata/Areas select Europe" > /tmp/preseed.txt
    echo "tzdata tzdata/Zones/Europe select Paris" >> /tmp/preseed.txt
    debconf-set-selections /tmp/preseed.txt

    apt-get update && apt-get install -y \
        openjdk-11-jre \
        wget \
        unzip \
        locales \
        perl \
        make \
        gcc \
        build-essential \
        liblocal-lib-perl \
        cpanminus

    # Configurer les locales
    locale-gen fr_FR.UTF-8
    update-locale LANG=fr_FR.UTF-8

    # Installer FastQC
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
    unzip fastqc_v0.11.9.zip
    chmod +x FastQC/fastqc
    mv FastQC /opt/FastQC

    # Installer le module Perl FindBin
    cpanm FindBin

%environment
    export PATH=/opt/FastQC:$PATH
    export LANG=fr_FR.UTF-8
    export LC_ALL=fr_FR.UTF-8
```

## `01_fastqc/sh.sh`

Source brute : lignes 48–107 de la compilation.

```bash
#!/bin/bash
if command -v sbatch &> /dev/null; then
 #SBATCH --job-name=fastqc
 #SBATCH --output=/dev/null
 #SBATCH --nodes=1
 #SBATCH --ntasks=48
 module load bioinfo/FastQC/0.12.1
 if true; then
  FILES=()
  for FILE in ../../RAW/*.fastq*; do
   FILENAME=$(basename "$FILE")
   FILENAME=$(echo "$FILENAME" | sed 's/_.*//')
   FILES+=("$FILENAME")
  done
  FILES=($(echo "${FILES[@]}" | tr ' ' '\n' | sort -u))
  n=6
  i=0
  j=0
  for ((i=0; i<${#FILES[@]}; i+=n)); do
   files=()
   for ((j=i; j<i+n && j<${#FILES[@]}; j++)); do
    files+=("${FILES[j]}")
   done
   for file in "${files[@]}"; do
    echo "${file}"
    (
     fastqc -t 4 -o ./ ../../RAW/${file}_1.fastq* &
     fastqc -t 4 -o ./ ../../RAW/${file}_2.fastq*
    ) &
   done
   wait
   chmod 777 *
  done
 fi
else
 if true; then
  FILES=()
  for FILE in ../../RAW/*.fastq*; do
   FILENAME=$(basename "$FILE")
   FILENAME=$(echo "$FILENAME" | sed 's/_.*//')
   FILES+=("$FILENAME")
  done
  FILES=($(echo "${FILES[@]}" | tr ' ' '\n' | sort -u))
  n=6
  i=0
  j=0
  for ((i=0; i<${#FILES[@]}; i+=n)); do
   files=()
   for ((j=i; j<i+n && j<${#FILES[@]}; j++)); do
    files+=("${FILES[j]}")
   done
   for file in "${files[@]}"; do
    echo "${file}"
    singularity exec ../fastqc.sif fastqc -t 4 -o ./ ../../RAW/${file}_1.fastq*
    singularity exec ../fastqc.sif fastqc -t 4 -o ./ ../../RAW/${file}_2.fastq*
   done
   chmod 777 *
  done
 fi
fi
```

## `02_multiqc/multiqc.def`

Source brute : lignes 113–140 de la compilation.

```text
Bootstrap: docker
From: ubuntu:20.04

%post
    export DEBIAN_FRONTEND=noninteractive

    # Installer les packages nécessaires
    apt-get update && apt-get install -y \
        python3 \
        python3-pip \
        wget \
        unzip \
        locales \
        perl \
        make \
        gcc \
        build-essential

    # Configurer les locales
    locale-gen fr_FR.UTF-8
    update-locale LANG=fr_FR.UTF-8

    # Installer MultiQC
    pip3 install multiqc

%environment
    export LANG=fr_FR.UTF-8
    export LC_ALL=fr_FR.UTF-8
```

## `02_multiqc/sh.sh`

Source brute : lignes 146–156 de la compilation.

```bash
#!/bin/bash
if command -v sbatch &> /dev/null; then
 #SBATCH --job-name=multiqc
 #SBATCH --outdir=/dev/null
 #SBATCH --nodes=1
 #SBATCH --ntasks=1
 module load bioinfo/MultiQC/1.14
 multiqc --outdir ./ ../*fastqc
else
 singularity exec ../multiqc.sif multiqc --outdir ./ ../*fastqc
fi
```

## `03_cutadapt/cutadapt.def`

Source brute : lignes 162–189 de la compilation.

```text
Bootstrap: docker
From: ubuntu:20.04

%post
    export DEBIAN_FRONTEND=noninteractive

    # Installer les packages nécessaires
    apt-get update && apt-get install -y \
        python3 \
        python3-pip \
        wget \
        unzip \
        locales \
        perl \
        make \
        gcc \
        build-essential

    # Configurer les locales
    locale-gen fr_FR.UTF-8
    update-locale LANG=fr_FR.UTF-8

    # Installer Cutadapt
    pip3 install cutadapt

%environment
    export LANG=fr_FR.UTF-8
    export LC_ALL=fr_FR.UTF-8
```

## `03_cutadapt/sh.sh`

Source brute : lignes 195–232 de la compilation.

```bash
#!/bin/bash
if command -v sbatch &> /dev/null; then
 #SBATCH --job-name=cutadapt
 #SBATCH --output=/dev/null
 #SBATCH --nodes=1
 #SBATCH --ntasks=1
 module load bioinfo/Cutadapt/4.3
 adapter1=$(cat ../../RAW/adapter1*.txt | tr -d '\n' | tr -d '\r')
 adapter2=$(cat ../../RAW/adapter2*.txt | tr -d '\n' | tr -d '\r')
 inputpath=/work/user/rlioutaud/RAW/hcon/rnaseq
 fastqs=$(ls $inputpath | grep -v '\.sh$' | sed 's/_[^_]*$//' | sort | uniq)
 for fastq in $fastqs; do
  samplename=$(basename "$fastq" | sed 's/_1\.fastq$//' | sed 's/_2\.fastq$//')
  cutadapt \
  -a $adapter1 \
  -A $adapter2 \
  -o ${samplename}_1.trim.fastq \
  -p ${samplename}_2.trim.fastq \
  $inputpath/${fastq}_1.fastq \
  $inputpath/${fastq}_2.fastq &
 done
 wait
else
 adapter1=$(cat ../../RAW/adapter1*.txt | tr -d '\n' | tr -d '\r')
 adapter2=$(cat ../../RAW/adapter2*.txt | tr -d '\n' | tr -d '\r')
 inputpath=/work/user/rlioutaud/RAW/hcon/rnaseq
 fastqs=$(ls $inputpath | grep -v '\.sh$' | sed 's/_[^_]*$//' | sort | uniq)
 for fastq in $fastqs; do
  samplename=$(basename "$fastq" | sed 's/_1\.fastq$//' | sed 's/_2\.fastq$//')
  singularity exec ../cutadapt.sif cutadapt \
  -a $adapter1 \
  -A $adapter2 \
  -o ${samplename}_1.trim.fastq \
  -p ${samplename}_2.trim.fastq \
  $inputpath/${fastq}_1.fastq \
  $inputpath/${fastq}_2.fastq
 done
fi
```

## `04_hisat2/autre/kallisto/sh.sh`

Source brute : lignes 240–318 de la compilation.

```bash
#!/bin/bash
# Chargement de kallisto
module load bioinfo/kallisto/0.50.1
# Emplacements utilisés
PathGenome=~/work/Hcontortus/rsrc/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa
PathRnaseq=../../0xRawData
PathTranscriptome=~/work/Hcontortus/rsrc/haemonchus_contortus.PRJEB506.WBPS18.mRNA_transcripts.fa
PathGff3=~/work/Hcontortus/rsrc/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3
# Début
if false; then # Création de l'index des transcrits
	kallisto index --index TRANSCRIPTS.ixd $PathTranscriptome
	chmod 777 TRANSCRIPTS.ixd
	wait
fi
if false; then # Application de kallisto pour chaque rna-seq
	declare -A lane_sample_list
	while read -r LANE SAMPLE_ID; do
		lane_sample_list["$LANE"]="$SAMPLE_ID"
	done < ./lane_sample_IDs.list
	for LANE in "${!lane_sample_list[@]}"; do
		SAMPLE_ID="${lane_sample_list[$LANE]}"
		kallisto quant \
             	 --index ./TRANSCRIPTS.ixd \
             	 --output-dir kallisto_${SAMPLE_ID}_out \
             	 --bootstrap-samples 100 \
             	 --threads 7 \
             	 ${PathRnaseq}/${LANE}_1*.fastq ${PathRnaseq}/${LANE}_2*.fastq >> kallisto_${SAMPLE_ID}.log 2>&1;
	done
	echo "Toutes les tâches kallisto sont terminées."
fi
if false; then # Déverrouillage de tous fichiers et suppression des log
	for elem in *; do
	    chmod -R 777 "$elem"
	done
	rm -r kallisto_SAMPLE_ID*
	#rm kallisto_*.log
	echo "Fichiers kallisto avec tous les droits."
fi
if false; then
	# Rangement des LOG
	mkdir kallisto_log && chmod 777 kallisto_log
	mv kallisto*.log ./kallisto_log/
	# Création des dossiers TSV H5 JSON
	mkdir kallisto_tsv && chmod 777 kallisto_tsv
	mkdir kallisto_h5 && chmod 777 kallisto_h5
	mkdir kallisto_json && chmod 777 kallisto_json
	for kallisto_folder in kallisto_*_out; do
		# Récupérer le nom du dossier sans extension
		kall_out_fold_name="${kallisto_folder%%_out}"
		# TSV
		# Renommages
                mv "${kallisto_folder}/abundance.tsv" "${kallisto_folder}/${kall_out_fold_name}_abundance.tsv"
                # Déplacement vers dossier
		mv "${kallisto_folder}/${kall_out_fold_name}_abundance.tsv" "kallisto_tsv/"
		# H5
		# Renommages
                mv "${kallisto_folder}/abundance.h5" "${kallisto_folder}/${kall_out_fold_name}_abundance.h5"
		# Déplacement vers dossier
                mv "${kallisto_folder}/${kall_out_fold_name}_abundance.h5" "kallisto_h5/"
		# JSON
		#Renommages
		mv "${kallisto_folder}/run_info.json" "${kallisto_folder}/${kall_out_fold_name}_run_info.json"
		# Déplacement vers dossier
		mv "${kallisto_folder}/${kall_out_fold_name}_run_info.json" "kallisto_json/"
        done
	rm -r kallisto_*_out
fi
if true; then # Récupération des données de chaque TSV puis concaténation en un fichier counts
	cd kallisto_tsv
	for tsv in *.tsv; do
		filename=$(basename -- "$tsv")
		filename="${filename%.*}"
		cut -f4 -d$'\t' $tsv | tail -n +2 > ${filename}.part.tsv
		echo -e "$filename" | cat - ${filename}.part.tsv > "${filename}.tmp.tsv"
		mv "${filename}.tmp.tsv" "${filename}.part.tsv"
	done
	paste *.part.tsv > counts.tsv
	cd ..
fi
```

## `04_hisat2/hisat2.def`

Source brute : lignes 324–356 de la compilation.

```text
Bootstrap: docker
From: ubuntu:20.04

%post
    export DEBIAN_FRONTEND=noninteractive

    # Installer les packages nécessaires
    apt-get update && apt-get install -y \
        python3 \
        python3-pip \
        wget \
        unzip \
        locales \
        perl \
        make \
        gcc \
        build-essential

    # Configurer les locales
    locale-gen fr_FR.UTF-8
    update-locale LANG=fr_FR.UTF-8

    # Installer HISAT2
    wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download -O hisat2.zip
    unzip hisat2.zip -d /opt/
    mv /opt/hisat2-2.2.0 /opt/hisat2

    # Ajouter HISAT2 au PATH
    ln -s /opt/hisat2/hisat2* /usr/local/bin/

%environment
    export LANG=fr_FR.UTF-8
    export LC_ALL=fr_FR.UTF-8
```

## `04_hisat2/sh.sh`

Source brute : lignes 362–519 de la compilation.

```bash
#!/bin/bash
# Path Genome
PathGenome=$(find "../../RAW/" -name "*.genomic.fa*" -print -quit)
Path=$(find ../../RES -type d -name "*cutadapt*" -print)
if command -v sbatch &> /dev/null; then
 #SBATCH --job-name=hisat2_alignment
 #SBATCH --output=/dev/null
 #SBATCH --nodes=1
 #SBATCH --ntasks=8
 module load devel/python/Python-3.11.1
 module load bioinfo/HISAT2/2.2.1
 module load bioinfo/samtools/1.19
 if true; then # Générer l'index
  IndexName="hisat2idx"
  hisat2-build -p 8 $PathGenome $IndexName
  chmod 777 hisat2idx*
 fi
 if true; then # Pour chaque fastq paired-end, génération d'un alignement SAM unique
  # Liste
  allfastqs=()
  for file in $Path/*.fastq; do
   filename=$(basename "$file")
   filename=$(echo "$filename" | sed 's/_.*//')
   allfastqs+=("$filename")
  done
  allfastqs=($(echo "${allfastqs[@]}" | tr ' ' '\n' | sort -u))
  echo "allfastqs ${allfastqs[@]}"
  # Sous-listes
  n=4
  i=0
  j=0
  for ((i=0; i<${#allfastqs[@]}; i+=n)); do
   fastqs=()
   for ((j=i; j<i+n && j<${#allfastqs[@]}; j++)); do
    fastqs+=("${allfastqs[j]}")
   done
   echo "fastqs ${fastqs[@]}"
   for fastq in "${fastqs[@]}"; do
    echo "fastq $fastq"
    if [[ "$fastq" == *""* ]]; then
     # Génère .SAM
     FirstFastq="${Path}/${fastq}_1.trim.fastq"
     SecondFastq="${Path}/${fastq}_2.trim.fastq"
     (
     hisat2 \
     -p 8 \
     --rg ID:"$fastq" \
     --rg SM:"$fastq" \
     --rg LB:"library" \
     --rg PL:"illumina" \
     -x hisat2idx \
     -1 ${FirstFastq} \
     -2 ${SecondFastq} \
     -S "./${fastq}.sam" &&
     # SAM -> BAM SORTED
     samtools view -bS "${fastq}.sam" | samtools sort -@ 4 -o "${fastq}.sorted.bam" &&
     # Supprime le SAM volumineux
     rm "${fastq}.sam"
     ) &
    fi
   done
   wait
  done
 fi
 if true; then
  for bam in *.sorted.bam; do
   samtools index -@ 4 "$bam"
  done
 fi
 if true; then # Contrôles qualité FlagStat
  echo "FlagStat ..."
  for BAM_FILE in *.sorted.bam; do
   BASENAME=$(basename "$BAM_FILE" .sorted.bam)
   samtools flagstat "$BAM_FILE" > "${BASENAME}_flagstat.txt"
  done
  echo "... terminé"
 fi
 if true; then # Contrôles qualité Stat
  echo "Stat ..."
  for BAM_FILE in *.sorted.bam; do
   BASENAME=$(basename "$BAM_FILE" .sorted.bam)
   samtools stats "$BAM_FILE" > "${BASENAME}_stats.txt"
  done
  echo "... terminé"
 fi
else
 if true; then # Générer l'index
  IndexName="hisat2idx"
  singularity exec ../hisat2.sif hisat2-build -p 8 $PathGenome $IndexName
  chmod 777 hisat2idx*
 fi
 if true; then # Pour chaque fastq paired-end, génération d'un alignement SAM unique
  # Liste
  allfastqs=()
  for file in $Path/*.fastq; do
   filename=$(basename "$file")
   filename=$(echo "$filename" | sed 's/_.*//')
   allfastqs+=("$filename")
  done
  allfastqs=($(echo "${allfastqs[@]}" | tr ' ' '\n' | sort -u))
  echo "allfastqs ${allfastqs[@]}"
  # Sous-listes
  n=4
  i=0
  j=0
  for ((i=0; i<${#allfastqs[@]}; i+=n)); do
   fastqs=()
   for ((j=i; j<i+n && j<${#allfastqs[@]}; j++)); do
    fastqs+=("${allfastqs[j]}")
   done
   echo "fastqs ${fastqs[@]}"
   for fastq in "${fastqs[@]}"; do
    echo "fastq $fastq"
    if [[ "$fastq" == *""* ]]; then
     # Génère .SAM
     FirstFastq="${Path}/${fastq}_1.trim.fastq"
     SecondFastq="${Path}/${fastq}_2.trim.fastq"
     singularity exec ../hisat2.sif hisat2 \
     -p 8 \
     --rg ID:"$fastq" \
     --rg SM:"$fastq" \
     --rg LB:"library" \
     --rg PL:"illumina" \
     -x hisat2idx \
     -1 ${FirstFastq} \
     -2 ${SecondFastq} \
     -S "./${fastq}.sam" &&
     # SAM -> BAM SORTED
     singularity exec ../samtools.sif samtools view -bS "${fastq}.sam" | singularity exec ../samtools.sif samtools sort -@ 4 -o "${fastq}.sorted.bam" &&
     # Supprime le SAM volumineux
     rm "${fastq}.sam"
    fi
   done
   wait
  done
 fi
 if true; then
  for bam in *.sorted.bam; do
   singularity exec ../samtools.sif samtools index -@ 4 "$bam"
  done
 fi
 if true; then # Contrôles qualité FlagStat
  echo "FlagStat ..."
  for BAM_FILE in *.sorted.bam; do
   BASENAME=$(basename "$BAM_FILE" .sorted.bam)
   singularity exec ../samtools.sif samtools flagstat "$BAM_FILE" > "${BASENAME}_flagstat.txt"
  done
  echo "... terminé"
 fi
 if true; then # Contrôles qualité Stat
  echo "Stat ..."
  for BAM_FILE in *.sorted.bam; do
   BASENAME=$(basename "$BAM_FILE" .sorted.bam)
   singularity exec ../samtools.sif samtools stats "$BAM_FILE" > "${BASENAME}_stats.txt"
  done
  echo "... terminé"
 fi
fi
```

## `05_MarkDuplicates/markduplicates.sh`

Source brute : lignes 525–558 de la compilation.

```bash
#!/bin/bash
#
# BAM triés
input=$(find ../../RES -type d -name "*hisat2*" -print)
#
echo "Marquage des duplicats ..."
module load bioinfo/Sambamba/1.0.1
# Liste des BAMs
allbams=()
for bam in $input/*.sorted.bam; do
 bamname=$(basename "$bam")
 bamname="${bamname%.*}"
 bamname="${bamname%.*}"
 allbams+=("$bamname")
done
allbams=($(echo "${allbams[@]}" | tr ' ' '\n' | sort -u))
# Mark duplicates
n=3
i=0
j=0
for ((i=0; i<${#allbams[@]}; i+=n)); do
 bams=()
 for ((j=i; j<i+n && j<${#allbams[@]}; j++)); do
  bams+=("${allbams[j]}")
 done
 for bam in "${bams[@]}"; do
  (
  sambamba markdup -t=8 "${input}/${bam}.sorted.bam" "./${bam}.mkdup.sorted.bam"
  echo "... one marked ..."
  ) &
 done
 wait
done
echo "... terminé"
```

## `05_MarkDuplicates/quality.sh`

Source brute : lignes 564–646 de la compilation.

```bash
#!/bin/bash
module load bioinfo/samtools/1.19
#
if false; then # Contrôles qualité FlagStat
 echo "FlagStat ..."
 for BAM_FILE in *.mkdup.sorted.bam; do
  (
   BASENAME=$(basename "$BAM_FILE" .mkdup.sorted.bam)
   samtools flagstat "$BAM_FILE" > "${BASENAME}_flagstat.txt"
  ) &
 done
 wait
 echo "... terminé"
fi
if false; then # Contrôles qualité Stat
 echo "Stat ..."
 for BAM_FILE in *.mkdup.sorted.bam; do
  (
   BASENAME=$(basename "$BAM_FILE" .mkdup.sorted.bam)
   samtools stats "$BAM_FILE" > "${BASENAME}_stats.txt"
  ) &
 done
 wait
 echo "... terminé"
fi
if true; then # Résumé qualité
 echo "Résumé qualité ..."
 # Création du fichier CSV avec les colonnes souhaitées
 output_csv="summary_metrics.csv"
 echo "File,Total Reads,Mapped Reads,Alignment Rate,Uniquely Mapped Reads,Duplication Rate,Coverage Depth" > "$output_csv"
 # Extraction et calcul des métriques
 for BAM_FILE in *.mkdup.sorted.bam; do
  BASENAME=$(basename "$BAM_FILE" .mkdup.sorted.bam)
  flagstat_file="${BASENAME}_flagstat.txt"
  stats_file="${BASENAME}_stats.txt"
  #
  # Extraction des données avec valeurs par défaut pour les champs vides
  total_reads=$(grep "in total" "$flagstat_file" | awk '{print $1}' | tr -d '\n\r')
  mapped_reads=$(grep "mapped (" "$flagstat_file" | head -n 1 | awk '{print $1}' | tr -d '\n\r')
  if [[ -n "$total_reads" && -n "$mapped_reads" ]]; then
    alignment_rate=$(echo "scale=4; $mapped_reads / $total_reads * 100" | bc)
  else
    alignment_rate=""
  fi
  secondary=$(grep "secondary" "$flagstat_file" | awk '{print $1}' | tr -d '\n\r')
  if [[ -n "$mapped_reads" && -n "$secondary" && -n "$total_reads" ]]; then
    uniquely_mapped_reads=$(echo "scale=4; (($mapped_reads - $secondary) / $total_reads) * 100" | bc)
  else
    uniquely_mapped_reads=""
  fi
  duplication_rate=$(grep "duplicates" "$flagstat_file" | head -n 1 | awk '{print $1}' | tr -d '\n\r')
  if [[ -n "$total_reads" && -n "$duplication_rate" ]]; then
    duplication_rate=$(echo "scale=4; $duplication_rate / $total_reads * 100" | bc)
  else
    duplication_rate=""
  fi
  rrna_rate=$(grep "QC-passed reads" "$flagstat_file" | awk '{print $1}' | tr -d '\n\r')
  if [[ -n "$total_reads" && -n "$rrna_rate" ]]; then
    rrna_rate=$(echo "scale=4; $rrna_rate / $total_reads" | bc)
  else
    rrna_rate=""
  fi
  genome_size=$(awk '{sum += $2} END {print sum}' ../../RAW/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa.fai)
  # Obtenir le nombre total de bases alignées à partir du fichier stats
  total_bases_aligned=$(grep "^SN" "$stats_file" | grep "bases mapped (cigar)" | head -n 1 | awk '{print $5}' | tr -d '\n\r')
  # Calculer la profondeur de couverture
  if [[ -n "$total_bases_aligned" && -n "$genome_size" ]]; then
   coverage_depth=$(echo "scale=4; $total_bases_aligned / $genome_size" | bc)
   coverage_depth="${coverage_depth}X"
  else
   coverage_depth=""
  fi
  exonic_rate=$(grep "QC-failed reads" "$flagstat_file" | awk '{print $1}' | tr -d '\n\r')
  if [[ -n "$total_reads" && -n "$exonic_rate" ]]; then
    exonic_rate=$(echo "scale=4; $exonic_rate / $total_reads" | bc)
  else
    exonic_rate=""
  fi
  # Écriture des résultats dans le fichier CSV sur une seule ligne
  echo "${BASENAME}.sorted.bam,${total_reads},${mapped_reads},${alignment_rate},${uniquely_mapped_reads},${duplication_rate},${coverage_depth}" >> "$output_csv"
 done
 echo "... terminé"
fi
```

## `05_MarkDuplicates/verifmarkduplicates.sh`

Source brute : lignes 652–657 de la compilation.

```bash
#!/bin/bash
module load bioinfo/samtools/1.19
#
file="T.R.ARA.F.1.mkdup.sorted.bam"
#
samtools flagstat $file
```

## `11_featurecounts/autres/htseq/sh.sh`

Source brute : lignes 666–669 de la compilation.

```bash
#!/bin/sh
if false; then # Compter les lectures alignées avec HTSeq
        #htseq-count -f bam -r pos -s no -i gene_id output.bam reference.gtf > counts.txt
fi
```

## `11_featurecounts/autres/mmquant/sh.sh`

Source brute : lignes 676–684 de la compilation.

```bash
#!/bin/bash
module load bioinfo/mmquant/1.0.9
#
annotation="../featureCounts/annotations.gtf"
for bam in ../../2xAlignement/hisat2/*.bam; do
	filename=$(basename -- "$bam")
	filename="${filename%.*}"
	mmquant --bam $bam --annotation $annotation --output ${filename}_counts.csv
done
```

## `11_featurecounts/autres/r.R`

Source brute : lignes 690–807 de la compilation.

```r
# Libraries
#library(DESeq2)
#library(tximport)
#library(GenomicFeatures)
#library(tidyverse)
library(ggplot2)
#
# Récupération des emplacements des fichiers d'abondance .tsv
if (!file.exists("synthesis.csv")) {
	data <- read.table("./rnaseq_metadata.txt", header = TRUE, sep = "\t")
	abundances <- data$path
	abundances <- as.vector(abundances)
	# Conteneur des résultats
	analyse <- list()
	# Préparation
	for (abundance in abundances) {
		# Récupération des données du fichier d'abondance
		theabundance <- read.table(abundance, header = TRUE, row.names = 1, sep = "\t")
		theabundance$gene <- rownames(theabundance)
		rownames(theabundance) <- NULL
		# Récupération du nom de fichier et ajout au tableau
		samplename <- sub(".*/([^/]+)/.*", "\\1", abundance)
		theabundance <- cbind(souche = samplename, theabundance)
		# Récupération du statut résistant et ajout au tableau
		resistant <- FALSE
		if (grepl("_ARA_", samplename) || grepl("_BET_", samplename) || grepl("_BUN_", samplename) || grepl("_MOU_", samplename)) {
			resistant <- TRUE
		}
		theabundance <- cbind(resistant = resistant, theabundance)
		# Concaténation du tableau au tableau général
		analyse <- rbind(analyse, theabundance)
	}
	df <- as.data.frame(analyse)
	# Liste des gènes
	genes <- unique(df$gene)
	#
	# Créer un tableau de 2 colonnes
	synthesis <- data.frame(
	  gene = c(),
	  se = c(),
	  re = c(),
	  del = c()
	)
	#
	for (thegene in genes) {
		subdf <- df[df$gene == thegene, ]
		sendf <- subdf[subdf$resistant == FALSE, ]
		resdf <- subdf[subdf$resistant == TRUE, ]
		sensible <- as.numeric(mean(sendf$tpm))
		resistant <- as.numeric(mean(resdf$tpm))
		delta <- as.numeric(resistant - sensible)
		newline <- data.frame(gene = thegene, se = sensible, re = resistant, del = delta)
		print(newline)
		synthesis <- rbind(synthesis, newline)
	}
	synthesis$del <- as.numeric(synthesis$del)
	synthesis <- synthesis[order(synthesis$del),]
	write.csv(synthesis, "./synthesis.csv", row.names = FALSE)
}
if (1 == 0) {
	library(readr)
	# Étape 1: Lire un fichier
	chemin_fichier <- "~/work/Hcontortus/rsrc/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3"
	donnees_brutes <- readLines(chemin_fichier)
	# Étape 2: Supprimer les lignes commençant par #
	donnees_filtrees <- donnees_brutes[!grepl("^#", donnees_brutes)]
	donnees_filtrees <- read.table(text = donnees_filtrees, sep = "\t", header = FALSE)
	donnees_filtrees <- donnees_filtrees[, 5:ncol(donnees_filtrees)]
	# Étape 3: Créer un dataframe
	donnees_df <- data.frame(accession = c(), gene = c())
	for (ligne in donnees_filtrees) {
		ac_partie_precedente <- gsub('.*\tID=transcript:(.*?);.*', '\\1', ligne)
		ac <- gsub(';.*', '', ac_partie_precedente)
		ge_partie_precedente <- gsub('.*description[=|:](.*?)[\\s|;].*', '\\1', ligne)
		ge <- gsub(';.*', '', ge_partie_precedente)
		if (length(ac) > 0 & length(ge) > 0) {
			newline <- data.frame(accession = ac, gene = ge)
			donnees_df <- rbind(donnees_df, newline)
			print(newline)
		}
	}
	# Étape 6: Exporter le dataframe en .csv
	write.csv(donnees_df, file = "/AccessionToGeneName.csv", row.names = FALSE)
}
if (file.exists("synthesis.csv")) {
	#if (0 == 1) {
	# Charger le fichier CSV en tant que dataframe
	synthesis <- read.csv("synthesis.csv", header = TRUE, sep = ",")
	# Charger la table de correspondances Annotation vs Fonction du gène
	AnnotationToGeneFunction <- as.data.frame(read.csv("~/work/Hcontortus/rsrc/AnnotationToGeneFunction.csv", header = FALSE, sep = ";"))
	# Boucle pour mettre à jour les valeurs de la colonne "gene" dans "synthesis"
	for (i in 1:nrow(synthesis)) {
	  gene_value <- synthesis$gene[i]
	  # Trouver l'index correspondant dans "Annotation" avec startsWith
	  index_annotation <- which(startsWith(gene_value, AnnotationToGeneFunction$X1))
	  # Vérifier si l'index est trouvé
	  if (length(index_annotation) > 0) {
	    # Mettre à jour la valeur dans "synthesis"
	    synthesis$gene[i] <- paste(gene_value, AnnotationToGeneFunction$X2[index_annotation], sep = " ")
	  }
	}
	# Calculer le nombre total de lignes
	nbrows <- nrow(synthesis)
	# Sélectionner les premiers et les derniers dixièmes des lignes
	top <- head(synthesis, nbrows / 250)
	end <- tail(synthesis, nbrows / 250)
	# Concaténer les deux parties
	synthesis <- rbind(top, end)
	#
	# Créer le graphique trié
        plot <- ggplot(synthesis, aes(x = reorder(gene, del), y = del)) +
          geom_bar(stat = "identity", fill = "blue") +
	  labs(x = "Gènes", y = "deltaTPM", title = "deltaTPM par gène") +
          theme(axis.text.x = element_text(size = 7, angle = 0, hjust = 1)) +
          coord_flip()
	# Sauvegarder le graphique trié en tant que fichier SVG
	ggsave("graph.svg", plot, device = "svg", height = 20)
}
```

## `11_featurecounts/autres/sh.sh`

Source brute : lignes 813–814 de la compilation.

```bash
#!/bin/sh
/work/user/rlioutaud/logiciels/R-4.3.2/bin/Rscript r.R
```

## `11_featurecounts/r.R`

Source brute : lignes 820–897 de la compilation.

```r
###########
# Libraries
#library(rtracklayer)
#library(GenomicFeatures)
library(GenomicRanges)
library(Rsubread)

#######
# Paths
annotations <- list.files("../../RAW", pattern = "\\.gff3$|\\.gff3\\.gz$|\\.gtf$|\\.gtf\\.gz$", full.names = TRUE)
gff3 <- annotations[grepl("\\.gff3$", annotations, ignore.case = TRUE)]
gtf <- annotations[grepl("\\.gtf$", annotations, ignore.case = TRUE)]
if (length(gtf) > 0) {
 gtf <- gtf[which.min(file.info(gtf)$size)]
}

###############
# featureCounts
if (FALSE) {# Traiter les BAMs
 # EMPLACEMENTS BAMs
 bams <- list.files("../../RES", pattern = ".*MarkDuplicates.*\\.bam$", full.names = TRUE, recursive = TRUE)
 # COMPTAGE
 counts <- featureCounts(
  files=bams,
  annot.ext="../../RAW/haemonchus_contortus.PRJEB506.WBPS18.canonical_geneset.gtf",
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  nthreads = 4
 )
 #Ecrire les résultats dans un fichier
 write.table(
  x=data.frame(counts$annotation[,c("GeneID","Length")],
  counts$counts,
  stringsAsFactors=FALSE),
  file="cnt.csv",
  quote=FALSE,
  sep=",",
  row.names=FALSE)
}
if (TRUE) {
 counts <- read.csv("cnt.csv", header=TRUE)
 sample_names <- colnames(counts)[-c(1, 2)]
 split_names <- strsplit(sample_names, "\\.")
 find_common_suffix <- function(split_names) {
  min_length <- min(sapply(split_names, length))
  common_suffix_length <- 0
  for (i in seq_len(min_length)) {
   suffix_parts <- sapply(split_names, function(x) x[length(x) - i + 1])
   if (length(unique(suffix_parts)) == 1) {
    common_suffix_length <- i
   } else {
    break
   }
  }
  return(common_suffix_length)
 }
 find_common_prefix <- function(split_names) {
  min_length <- min(sapply(split_names, length))
  common_prefix_length <- 0
  for (i in seq_len(min_length)) {
   prefix_parts <- sapply(split_names, function(x) x[i])
   if (length(unique(prefix_parts)) == 1) {
    common_prefix_length <- i
   } else {
    break
   }
  }
  return(common_prefix_length)
 }
 common_suffix_length <- find_common_suffix(split_names)
 common_prefix_length <- find_common_prefix(split_names)
 new_sample_names <- sapply(split_names, function(x) {
  paste(x[(common_prefix_length + 1):(length(x) - common_suffix_length)], collapse = ".")
 })
 colnames(counts)[-c(1, 2)] <- new_sample_names
 print(colnames(counts))
 write.csv(counts, "cnt.csv", row.names=FALSE, quote=FALSE)
}
```

## `11_featurecounts/sh.sh`

Source brute : lignes 903–905 de la compilation.

```bash
#!/bin/sh
/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript r.R
chmod 777 *
```

## `13_normalisation_cpm/cpm.r.sh`

Source brute : lignes 912–944 de la compilation.

```r
#!/bin/sh
#####################
# R executable path #
#####################
Rscript="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
##########
# Script #
##########
read -r -d '' scriptR <<'EOT'
# vvvvvvvv
# CPM fonction
library(tools)
cpm <- function(cnt) {
 filename <- tools::file_path_sans_ext(cnt)
 all <- read.csv(cnt, header = TRUE, row.names = 1)
 GeneID <- rownames(all)
 Length <- all[, 1]
 cnt <- as.matrix(all[, -1])
 TotalMapped <- colSums(cnt)
 CPM <- t(t(cnt) / TotalMapped) * 1e6
 CPM <- cbind(GeneID, Length, CPM)
 CPM <- CPM[order(rownames(CPM)), ]
 write.csv(CPM, file = paste(filename, ".cpm.csv", sep=""), row.names = FALSE, quote = FALSE)
}
for (file in file.path(list.files("../", pattern = ".*featurecounts.*", full.names = TRUE), "cnt.csv")) {
 cpm(file)
}
# ^^^^^^^^
EOT
#####################
# Rscript exécution #
#####################
echo "$scriptR" | $Rscript -
```

## `13_normalisation_mrm/mrm.r.sh`

Source brute : lignes 950–994 de la compilation.

```r
#!/bin/sh
#####################
# R executable path #
#####################
Rscript="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
##########
# Script #
##########
read -r -d '' scriptR <<'EOT'
# vvvvvvvv
# MPM Mean of Ratio Method de DESeq2 fonction
library(tools)
library(DESeq2)
mrm <- function(cnt) {
 filename <- tools::file_path_sans_ext(cnt)
 all <- read.csv(cnt, header = TRUE, row.names = 1)
 GeneID <- rownames(all)
 Length <- all[, 1]
 cnt <- as.matrix(all[, -1])
 #
 rownames(cnt) <- GeneID
 colData <- data.frame(condition = colnames(cnt))
 design <- ~ 1
 #
 dds <- DESeqDataSetFromMatrix(countData = cnt, colData = colData, design = design)
 # Normalisation interne de DESeq2
 dds <- estimateSizeFactors(dds)
 # Extraction des comptages normalisés
 MRM <- counts(dds, normalized = TRUE)
 # Export des comptages normalisés
 MRM <- cbind(GeneID, Length, MRM)
 MRM <- MRM[order(rownames(MRM)), ]
 write.csv(MRM, file = paste("./", basename(filename), ".mrm.csv", sep = ""), row.names = FALSE, quote = FALSE)
}

for (file in file.path(list.files("../", pattern = ".*featurecounts.*", full.names = TRUE), "cnt.csv")) {
 mrm(file)
}
# ^^^^^^^^
EOT
#####################
# Rscript exécution #
#####################
echo "$scriptR" | $Rscript -
chmod 777 *
```

## `13_normalisation_tmm/tmm.r.sh`

Source brute : lignes 1001–1038 de la compilation.

```r
#!/bin/sh
#####################
# R executable path #
#####################
Rscript="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
##########
# Script #
##########
read -r -d '' scriptR <<'EOT'
# vvvvvvvv
# TMM déCPMisé fonction
library(tools)
library(edgeR)
tmm <- function(cnt) {
 filename <- tools::file_path_sans_ext(cnt)
 all <- read.csv(cnt, header = TRUE, row.names = 1)
 GeneID <- rownames(all)
 Length <- all[, 1]
 cnt <- as.matrix(all[, -1])
 TotalMapped <- colSums(cnt)
 dge <- DGEList(counts = cnt)
 dge <- calcNormFactors(dge, method = "TMM")
 TMM <- cpm(dge, normalized.lib.sizes = TRUE)
 colnames(TMM) <- colnames(cnt)
 TMM <- TMM * TotalMapped / 1e6
 TMM <- cbind(GeneID, Length, TMM)
 TMM <- TMM[order(rownames(TMM)), ]
 write.csv(TMM, file = paste(filename, ".tmm.csv", sep=""), row.names = FALSE, quote = FALSE)
}
for (file in file.path(list.files("../", pattern = ".*featurecounts.*", full.names = TRUE), "cnt.csv")) {
 tmm(file)
}
# ^^^^^^^^
EOT
#####################
# Rscript exécution #
#####################
echo "$scriptR" | $Rscript -
```

## `13_normalisation_tpm/tpm.r.sh`

Source brute : lignes 1044–1083 de la compilation.

```r
#!/bin/sh
#####################
# R executable path #
#####################
Rscript="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
##########
# Script #
##########
read -r -d '' scriptR <<'EOT'
# vvvvvvvv
# TPM fonction
library(tools)
tpm <- function(cnt) {
 filename <- tools::file_path_sans_ext(cnt)
 all <- read.csv(cnt, header = TRUE, row.names = 1)
 GeneID <- rownames(all)
 Length <- all[, 1]
 cnt <- as.matrix(all[, -1])
 for (col in colnames(cnt)) {
  for (lin in 1:nrow(cnt)) {
   cnt[lin, col] <- cnt[lin, col] / Length[lin]
  }
  colsum <- sum(cnt[, col])
  for (lin in 1:nrow(cnt)) {
   cnt[lin, col] <- cnt[lin, col] / colsum * 1e6
  }
 }
 tpms <- cbind(GeneID, Length, cnt)
 tpms <- tpms[order(rownames(tpms)), ]
 write.csv(tpms, file = paste("./", basename(filename), ".tpm.csv", sep=""), row.names = FALSE, quote = FALSE)
}
for (file in file.path(list.files("../", pattern = ".*featurecounts.*", full.names = TRUE), "cnt.csv")) {
 tpm(file)
}
# ^^^^^^^^
EOT
#####################
# Rscript exécution #
#####################
echo "$scriptR" | $Rscript -
```

## `14_normalisation/brouillon_autres_normalisations.r`

Source brute : lignes 1090–1151 de la compilation.

```r
library(DGEobj.utils)
path <- "../3vCounting/featureCounts/counts.txt"
#####
# CNT
# Récupération de counts.txt de featureCounts et création de cnt.csv
GeneLengthsAndCounts <- read.delim(path)
#
GeneID <- GeneLengthsAndCounts[, 1]
Length <- GeneLengthsAndCounts[, 2]
NbOfReadsMappedToGenes <- GeneLengthsAndCounts[, -c(1, 2)]
TotalNbOfMappedReads <- colSums(NbOfReadsMappedToGenes)
cnt <- cbind(GeneID, Length, NbOfReadsMappedToGenes)
write.csv(cnt, file = "cnt.csv", row.names = FALSE)
########################
# TMM : faire ou non ? #
TMMcheck <- TRUE
########################

######
# RPKM
# fonction
RPKM <- function(CPM) {
 GeneID <- rownames(CPM)
 Length <- CPM[, 1]
 RPKMs <- CPM[, -1]
 for (colonne in names(RPKMs)) {
  colsum <- sum(RPKMs[[colonne]])
  for (ligne in 1:nrow(RPKMs)) {
   RPKMs[ligne, colonne] <- RPKMs[ligne, colonne] * (1000 / Length[ligne])
  }
 }
 #rpkm <- RPKMs * (1e3 / Length)
 rpkm <- cbind(GeneID, Length, RPKMs)
 return(rpkm)
}
# calcul
if (TMMcheck == FALSE) {
 cpm <- read.csv("cpm.csv", header = TRUE, row.names = 1)
} else {
 cpm <- read.csv("tmm_cpm.csv", header = TRUE, row.names = 1)
}
if (TRUE) {# RPKM, avec équation
 rpkm <- RPKM(cpm)
 if (TMMcheck == FALSE) {
  write.csv(rpkm, file = "rpkm.csv", row.names = FALSE)
 } else {
  write.csv(rpkm, file = "tmm_rpkm.csv", row.names = FALSE)
 }
}
if (FALSE) {
 dge <- DGEList(counts = NbOfReadsMappedToGenes)
 dge <- calcNormFactors(dge)
 normCounts <- cpm(dge)
 normCounts <- normCounts / TotalNbOfMappedReads
 for (i in 1:length(GeneLengths)) {
   normCounts[,i] <- normCounts[,i] / GeneLengths[i]
 }
 rpkm <- normCounts * 1e3
 TotalNbOfMappedReads <- colSums(rpkm)
 rpkm <- rpkm / TotalNbOfMappedReads * 1e6
 write.csv(rpkm, file = "rpkm_edgeR.csv", row.names = TRUE)
}
```

## `14_normalisation/sh.sh`

Source brute : lignes 1156–1161 de la compilation.

```bash
#!/bin/sh
Rpath="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
$Rpath tpm.r
#$Rpath tmm.r
#$Rpath cpm.r
#$Rpath mrm.r
```

## `15_exploration/sh.r.sh`

Source brute : lignes 1168–1422 de la compilation.

```r
#!/bin/sh
#####################
# R executable path #
#####################
Rscript="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
##########
# Script #
##########
read -r -d '' scriptR <<'EOT'
# vvvvvvvv
PCA <- function(cntfile, tf = NULL, getconditions = NULL) {
 library(tidyverse)
 library(plotly)
 library(htmltools)
 # Importation du fichier CSV et préparation des données
 df <- read.csv(cntfile, row.names = 1, header = TRUE)
 filename <- basename(cntfile)
 ExceptLastExt <- substr(
                   filename,
                   gregexpr("\\.", filename)[[1]][1] + 1,
                   gregexpr("\\.", filename)[[1]][length(gregexpr("\\.", filename)[[1]])] - 1
                  )
 #df <- head(df, 600)
 df <- df[, -1]
 #
 ##########################
 ##########################
 # Exclusion d'échantillons
 if (FALSE) {
  substrings_to_exclude <- c("nothing", "nothingelse")
  pattern <- paste(substrings_to_exclude, collapse = "|")
  cols_to_keep <- !sapply(names(df), function(x) grepl(pattern, x))
  df <- df[, cols_to_keep]
 }
 ##########################
 ##########################
 #
 df <- t(as.matrix(df))
 samples <- rownames(df)
 # Ajouter une petite constante aux données pour éviter les zéros
 epsilon <- 1# +1 aux counts, avant log2 transform, pcq log2(0)=infini
 df <- df + epsilon
 # Post-normalisation transformation
 if (!is.null(tf)) {
  if (tf == "rlog") {
   df <- log(df)
  } else if (tf == "vst") {
   variances <- apply(df, 2, var)
   vst_transform <- sqrt(variances)
   df <- log(df + 1) / vst_transform
  }
 }
 # Supprime colonnes (maintenant gènes) vides
 df <- df[, apply(df, 2, function(x) !all(is.na(x)) && length(unique(x)) > 1)]
 ########################################################
 # Effectuer une analyse en composantes principales (PCA)
 pca <- prcomp(df, scale. = TRUE)
 ########################################################
 #
 ###############################################
 # Récupération des conditions en liste de liste
 split_files <- strsplit(samples, "\\.")
 max_length <- max(sapply(split_files, length))
 conditions <- vector("list", max_length)
 for (i in seq_len(max_length)) {
   elements_at_i <- sapply(split_files, function(x) ifelse(length(x) >= i, x[i], NA))
   conditions[[i]] <- sort(unique(elements_at_i))
 }
 NbConditions <- length(conditions)
 print(conditions)
 # Récupérer le deuxième élément du troisième item
 second_element_third_item <- conditions[[3]][2]
 ###############################################
 #
 #
 ##################
 # Stock des plotly
 plots <- list()
 ##################
 # Variance expliquée
 plot_eigen <- plot_ly(x = 1:length(pca$sdev^2), y = pca$sdev^2, type = "bar") %>%
  layout(xaxis = list(title = "Composantes principales"),
         yaxis = list(title = "Variance expliquée"),
         title = "Valeurs propres des PCs"
        )
 plots[[length(plots) + 1]] <- plot_eigen
 # PCs
 NbPCs <- min(5, pca$rank)
 # % variance expliquée
 varexp <- summary(pca)$importance[2, ]
 varexp <- varexp^2
 totvar <- sum(varexp)
 varexp <- round(varexp / totvar * 100, digits=2)
 #
 # Pour chaque paire de PCs
 for (p in 1:NbPCs) {
  # Pour chaque condition
  for (condition in conditions) {
   # Récupération du nom de la condition
   cond=""
   for (case in condition) {
    cond <- paste0(cond, case)
   }
   #
   plot_pca <- plot_ly() %>%# on initialise un plotly,
    layout(
     xaxis = list(
      title = list(
       text = paste0("PC", p, " var ", varexp[p], " %"),
       font = list(size = 18, color = "black", weight = "bold"),
       standoff = 10
       )
      ),
     yaxis = list(
      title = list(
       text = paste0("PC", p+1, " var ", varexp[p+1], " %"),
       font = list(size = 18, color = "black", weight = "bold")
       )
      ),
     legend = list(
      orientation = 'h',
      y = 1.18,
      x = 0.5,
      xanchor = "center",
      font = list(size = 18, color = "black", weight = "bold")
      )
     )
   ExistingColors <- c("red", "blue", "green", "orange", "purple", "yellow", "cyan", "magenta", "black", "gray", "white")
   #
   # Pour chaque cas de la condition en cours de traitement
   i <- 1
   for (case in condition) {
    case <- case
    color <- ExistingColors[i]
    if (cond == 'FM') {
     if (case == 'F') {
      color <- "coral"
     }
     if (case == 'M') {
      color <- "darkcyan"
     }
    }
    if (cond == 'RS') {
     if (case == 'R') {
      color <- "red"
     }
     if (case == 'S') {
      color <- "green"
     }
    }
    pattern <- paste0("(^", case, "\\.)|(\\.", case, "$)|(\\.", case, "\\.)")
    pcasub <- pca$x[grepl(pattern, rownames(pca$x)), ]
    name <- case
    name <- gsub("\\bR\\b", "Resistant", name)
    name <- gsub("\\bS\\b", "Susceptible", name)
    name <- gsub("\\bM\\b", "Male", name)
    name <- gsub("\\bF\\b", "Female", name)
    plot_pca <- plot_pca %>% add_trace(
     x = pcasub[, paste0("PC", p)],
     y = pcasub[, paste0("PC", p+1)],
     mode = "markers",
     marker = list(color = color),
     text = rownames(pcasub),
     type = "scatter",
     textposition = "top",
     name = name)
    i <- i + 1
   }
   plot_pca <- plot_pca %>% config(toImageButtonOptions = list(format = 'svg', filename = 'pca'))
   plots[[length(plots) + 1]] <- plot_pca
  }
 }
 # Multiplot
 html_content <- tagList(div(style = "display: flex; flex-direction: column;",
                             lapply(plots, function(plot) div(plot))))
 # Sauvegarder le fichier HTML
 save_html(html_content, file = "plots.html")
 #
 library(pheatmap)
 library(RColorBrewer)
 # CONTROLE QUALITE CORRELATION SAMPLES HEATMAP
 sampleCor <- cor(t(df))
 sampleDists <- as.dist(1 - cor(t(df)))
 sampleDistMatrix <- as.matrix(sampleDists)
 # Echelle de la heatmap
 valmin <- min(sampleDistMatrix)
 valmax <- max(max(sampleDistMatrix), 0)
 ######################################################
 # il faut AUTANT de props que de couleurs en palette !
 echelle <- "brewer"
 #
 if (echelle == "kelvin") {
  palette <- c("darkblue", "cyan", "white", "yellow", "orange", "red")
  props <- seq(0, 1, length.out = length(palette))
 }
 if (echelle == "rainbow") {
  palette <- rev(c("purple", "darkblue", "cyan", "green", "yellow", "orange", "red"))
  props <- seq(0, 1, length.out = length(palette))
 }
 if (echelle == "thermic") {
  palette <- c("black", "purple", "magenta", "red", "orange", "yellow", "white")
  props <- seq(0, 1, length.out = length(palette))
 }
 if (echelle == "tricolore") {
  palette <- c("blue", "yellow", "red")
  props <- c(0, 0.5, 1)
 }
 if (echelle == "rvb") {
  palette <- rev(c("blue", "green", "red"))
  props <- c(0, 0.5, 1)
 }
 if (echelle == "brewer") {
  palette <- RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")
  #palette <- palette[1:5]
  props <- seq(0, 1, length.out = length(palette))
 }
 if (echelle == "viridis") {
  palette <- c("purple", "darkblue", "darkcyan", "mediumseagreen", "yellow")
  props <- seq(0, 1, length.out = length(palette))
 }
 ######################################################
 allcolors <- c()
 for (i in 1:(length(palette) - 1)) {
  num_colors_segment <- round((props[i + 1] - props[i]) * 100)
  palette_segment <- colorRampPalette(c(palette[i], palette[i + 1]))(num_colors_segment)
  allcolors <- c(allcolors, palette_segment)
 }
 #
 breaks <- valmin + props * (valmax - valmin)
 legend_breaks <- breaks
 #########
 # Heatmap
 svg("ClusteredHeatmap.svg")
 pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  border_color = NA,
  color = allcolors,
  breaks = seq(from = valmin, to = valmax, length.out = length(allcolors) + 1),
  legend_breaks = legend_breaks,
  legend_labels = round(legend_breaks, 3)
 )
 #color = allcolors,
 dev.off()
}
for (file in file.path(list.files("../", pattern = ".*normalisation_mrm.*", full.names = TRUE), ".*mrm.csv")) {
 PCA(file, tf="rlog")
}
# ^^^^^^^^
EOT
#####################
# Rscript exécution #
#####################
echo "$scriptR" | $Rscript -
```

## `16_deseq2/autres/edgeR/r.R`

Source brute : lignes 1431–1495 de la compilation.

```r
######################
# Chargement de DESeq2
library(edgeR)
library(ggplot2)
###########################
# Création d'un fichier log
###########################
log <- file("log.log", open = "w")
##############################################################################
writeLines("Chargement de la table de comptage donnée par featureCounts", log)
##############################################################################
counts <- read.table(
 "/save/user/rlioutaud/pipelines/rnaseq/3vCounting/featureCounts/counts.txt",
 header=TRUE,
 row.names=1,
 sep="\t"
)
counts <- counts[, -1]
#
writeLines(capture.output(head(counts)), log)
#####################################################################################################
writeLines("Récup° des conditions exp. et écriture dans table échantillon-conditions_multiples", log)
#####################################################################################################
SamplesNames <- colnames(counts)
SamplesConditionsAll <- list()
for (i in seq_along(SamplesNames)) {
 SampleConditionsAll <- as.vector(strsplit(SamplesNames[i], "\\.")[[1]])
 SamplesConditionsAll <- append(SamplesConditionsAll, list(SampleConditionsAll))
}
BorneGaucheConditions <- 2
#BorneDroiteConditions <- length(SamplesConditionsAll[[1]]) - 3
BorneDroiteConditions <- 4
coldata <- data.frame(name = SamplesNames)
# Pour autant qu'il y a de conditions
for (i in BorneGaucheConditions:BorneDroiteConditions) {
 ConditionValues <- c()
 # Pour chaque échantillon
 for (j in seq_along(SamplesNames)) {
  NewValue <- SamplesConditionsAll[[j]][i]
  ConditionValues <- c(ConditionValues, NewValue)
 }
 ConditionCategories <- sort(unique(ConditionValues))
 ConditionName <- paste(ConditionCategories, collapse = "")
 coldata[[ConditionName]] <- as.vector(factor(ConditionValues))
}
coldata[["type"]] <- rep("paired-end", nrow(coldata))
coldata$RS <- as.factor(coldata$RS)
writeLines(capture.output(head(coldata)), log)
##############################################
writeLines("Création d'un objet DGEList", log)
##############################################
dge <- DGEList(counts = counts, group = coldata$RS) # Condition principale
# Filtrer les données si nécessaire
# dge <- dge[rowSums(dge$counts) > 10, ]
# Normalisation (par exemple, TMM normalization)
dge <- calcNormFactors(dge)
# Créer un modèle de design avec des interactions entre les conditions
design <- model.matrix(~0 + RS + FM + RS:FM, data=coldata)
# Effectuer l'analyse différentielle en utilisant edgeR
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef=2)  # coefficient associé à l'effet de la condition croisée
# Récupérer les résultats
results <- topTags(qlf)
# Afficher les résultats
writeLines(capture.output(results), log)
```

## `16_deseq2/autres/edgeR/sh.sh`

Source brute : lignes 1501–1502 de la compilation.

```bash
#!/bin/sh
/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript r.R
```

## `16_deseq2/autres/sleuth/r.R`

Source brute : lignes 1509–1620 de la compilation.

```r
library(sleuth)
library(tidyverse)
library(patchwork)

# load metadata for the RNAseq reads processed by kallisto
metadata <- read.table("rnaseq_metadata.txt", header=T)

so <- sleuth_prep(metadata, extra_bootstrap_summary = TRUE, num_cores=2)

######
# PCAs

pca_matrix <- sleuth:::spread_abundance_by(so$obs_norm, "est_counts", so$sample_to_covariates$sample)

sample_pca <- prcomp(t(pca_matrix))

# "sdev" contains the standard deviation explained by each PC, so if we square it we get the eigenvalues (or explained variance)
# "rotation" contains the variable loadings for each PC, which define the eigenvectors
# "x" contains the PC scores, i.e. the data projected on the new PC axis
# "center" in this case contains the mean of each gene, which was subtracted from each value
# "scale" contains the value FALSE because we did not scale the data by the standard deviation 

pc_eigenvalues <- sample_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues




# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- sample_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

# print the result
pc_scores


pc_scores <- dplyr::left_join(pc_scores, so$sample_to_covariates, by = "sample")


plot_pca_pop <-
    ggplot(pc_scores, aes(x = PC1, y = PC2, col=population)) +
    geom_point() +
    theme_bw() +
    labs(x = paste0("PC1 variance: ",round(pc_eigenvalues$pct[1],digits=2),"%"),
        y = paste0("PC2 variance: ",round(pc_eigenvalues$pct[2],digits=2),"%"))

plot_pca_sex <-
    ggplot(pc_scores, aes(x = PC1, y = PC2, col=sex)) +
    geom_point() +
    theme_bw() +
    labs(x = paste0("PC1 variance: ",round(pc_eigenvalues$pct[1],digits=2),"%"),
        y = paste0("PC2 variance: ",round(pc_eigenvalues$pct[2],digits=2),"%"))

plot_pca_drug <-
    ggplot(pc_scores, aes(x = PC1, y = PC2, col=drug)) +
    geom_point() +
    theme_bw() +
    labs(x = paste0("PC1 variance: ",round(pc_eigenvalues$pct[1],digits=2),"%"),
        y = paste0("PC2 variance: ",round(pc_eigenvalues$pct[2],digits=2),"%"))


plot_pca_pop + plot_pca_sex + plot_pca_drug + plot_layout(ncol=3)

ggsave("figure_rnaseq_pca_pop_sex_drug.png")

###########
# 2e figure

abund <- sleuth:::spread_abundance_by(so$obs_norm, "tpm", so$sample_to_covariates$sample)

all_pairs <- sleuth:::apply_all_pairs(abund, sleuth:::jsd)

s2c <- so$sample_to_covariates

rownames(s2c) <- s2c$sample
annotation_cols = setdiff(colnames(so$sample_to_covariates),"sample")
s2c <- s2c[, annotation_cols, drop = FALSE]

color_high = "white"
color_low = "cornflowerblue"

    colors <- colorRampPalette(c(color_high, color_low))(100)
    
pdf(file = "figure_rnaseq_heatmap_pop_sex_drug.pdf")
    
pheatmap::pheatmap(all_pairs, annotation_col = s2c,
    color = colors, cluster_rows = TRUE, cluster_cols = TRUE,
    clustering_distance_cols = dist(all_pairs), clustering_distance_rows = dist(all_pairs),
    treeheight_row = 0)

 dev.off()

png(file = "figure_rnaseq_heatmap_pop_sex_drug.png")
    
pheatmap::pheatmap(all_pairs, annotation_col = s2c,
    color = colors, cluster_rows = TRUE, cluster_cols = TRUE,
    clustering_distance_cols = dist(all_pairs), clustering_distance_rows = dist(all_pairs),
    treeheight_row = 0)

 dev.off()
```

## `16_deseq2/autres/sleuth/sh.sh`

Source brute : lignes 1626–1627 de la compilation.

```bash
#!/bin/sh
/work/user/rlioutaud/logiciels/R-4.3.2/bin/Rscript r.R
```

## `16_deseq2/deseq2.r.sh`

Source brute : lignes 1633–1767 de la compilation.

```r
#!/bin/sh
#####################
# R executable path #
#####################
Rscript="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
##########
# Script #
##########
read -r -d '' scriptR <<'EOT'
# vvvvvvvv
# Chargement de DESeq2 #####
library(DESeq2)
library(ggplot2)
library(tools)
# Création d'un fichier log si n'existe pas déjà #####
if (!file.exists("log.log")) {
 log <- file("log.log", open = "w")
} else {
 log <- file("log.log", open = "a+")
}
logw <- function(textline) { writeLines(textline, log) }
# Sélection de la / des table de comptage par nom
path <- "../11_featurecounts/"
files <- list.files(path, full.names = TRUE)
files <- files[grep("cnt.csv", files)]
#
# Dans le dossier désigné, pour chaque fichier en .csv
for (file in files) {
 filename <- tools::file_path_sans_ext(basename(file))
 if (TRUE) { logw("Récupération des counts");
  writeLines("Chargement de la table de comptage donnée par featureCounts", log) #####
  counts <- read.csv(file, header=TRUE, row.names=1)
  counts <- counts[, -1]# Supprime la colonne Longueur des gènes
  counts <- as.matrix(counts)
  counts <- matrix(as.integer(counts), nrow = nrow(counts), ncol = ncol(counts), dimnames = list(rownames(counts), colnames(counts)))
  # Fonctions pour sélectionner et exclure des colonnes
  SeleCols <- function(data, pattern) {
   if (pattern != "") {
    cols <- grep(pattern, colnames(data))
   } else {
    cols <- 1:ncol(data)
   }
   return(data[, cols])
  }
  ExCols <- function(data, pattern) {
   if (pattern != "") {
    exclude_cols <- grep(pattern, colnames(data))
    data <- data[, -exclude_cols]
   }
   return(data)
  }
  # Filtrage de la table globale, "" pour ne pas garder/extraire, | pour séparer motifs multiples, \\. pour .
  counts <- ExCols(counts, "")
  counts <- SeleCols(counts, "\\.F\\.")
  print("counts")
  print(head(counts))
  # Création de la table de référence et filtrage des colonnes, "" pour ne pas garder/extraire
  reference <- counts
  reference <- ExCols(reference, "")
  reference <- SeleCols(reference, "\\.CHI\\.|\\.LUC\\.")
  print("reference")
  print(head(reference))
  # Création de la table de référence et filtrage des colonnes, "" pour ne pas garder/extraire
  comparison <- counts
  comparison <- ExCols(comparison, "")
  comparison <- SeleCols(comparison, "\\.BET\\.")
  print(head("comparison"))
  print(head(comparison))
  # Créer un DataFrame combiné
  CountsRefComp <- cbind(reference, comparison)
  # Créer un vecteur de conditions correspondant aux colonnes
  conditions <- factor(c(rep("ref", ncol(reference)), rep("comp", ncol(comparison))), levels=c("ref", "comp"))
  print(levels(conditions))
  # Créer un objet DESeq2DataSet
  dds <- DESeqDataSetFromMatrix(countData = CountsRefComp, colData = data.frame(condition = conditions), design = ~ condition)
  if (TRUE) { logw("Pre-filtering");
   smallestGroupSize <- 3
   keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
   dds <- dds[keep,]
  }
  ####################
  logw("Normalisation DESeq2")
  dds <- estimateSizeFactors(dds)
  if (TRUE) { logw("Export de la table normalisée vers le dossier source")
   normcnt <- counts(dds, normalized=TRUE)
   write.csv(normcnt, file = paste0(path, filename, ".mrm.csv"), row.names = TRUE)
  }
  logw("Analyse DGE DESeq2")
  # DESeq() détecte automatiquement la normalisation préalablement faite
  if (TRUE) { logw("DESeq() par défaut avec outlier replacement minRep");
   dds <- DESeq(dds)
  } else { logw("DESeq() spécifié sans outlier replacement minRep=Inf")
   # sans outlier = équivalent de estimateSizeFactors+estimateDispersions+nbinomWaldTest
   dds <- DESeq(dds, minRep=Inf)
  }
 }
 #########################  ^ Objet DDS  v Extraction & DGE  #########################
 logw("Extract° de l'objet DESeq2 d'un tab. des résult. puis Plots DGE qualité")
 res <- results(dds)
 res <- res[order(res$padj), ]
 write.csv(res, file = paste0(filename, ".deseq2.dge.csv"), row.names = TRUE)
 logw(capture.output(head(dds)))
 logw(capture.output(head(res)))
 logw("Génération des plot conventionnels en RNAseq DGE") #####
 if (TRUE) { logw("Plot de dispersion");
  pdf(paste0(filename, "_", "dispersion", ".pdf"))
  plotDispEsts(dds)
  dev.off()
 }
 if (TRUE) {
  svg(paste0(filename, "_", "dispersion", ".svg"))
  dispersions <- dispersionFunction(dds)
  #dispersions_df <- data.frame(dispersion = dispersions, group = group)
  #mean_dispersion <- aggregate(dispersion ~ group, dispersions_df, mean)$dispersion
  plot(dispersions, type = "b", col = "blue", pch = 16, cex = 1.5, main = "Custom Dispersion Plot")
  #points(1:length(mean_dispersion), mean_dispersion, type = "b", col = "red", pch = 16, lwd = 2)
  dev.off()
 }
 if (FALSE) { logw("Plot MA des résultats de l'analyse différentielle");
  pdf(paste0(filename, "_", "plot_MA.pdf"))
  plotMA(res, ylim=c(-2,2))
  dev.off()
 }
 # Suppression des MRM générés
 mrm <- list.files(path, full.names = TRUE)
 mrm <- mrm[grep(".mrm.csv", mrm)]
 file.remove(mrm)
 # Fin boucle for VVV
}
# ^^^^^^^^
EOT
#####################
# Rscript exécution #
#####################
echo "$scriptR" | $Rscript -
```

## `17_DGEplot/dgeplot.r.sh`

Source brute : lignes 1773–1949 de la compilation.

```r
#!/bin/sh
#####################
# R executable path #
#####################
Rscript="/save/user/rlioutaud/bin/R-4.3.2/bin/Rscript"
##########
# Script #
##########
read -r -d '' scriptR <<'EOT'
# vvvvvvvv
# Chargement de DESeq2 #####
library(ggplot2)
library(tools)
# Création d'un fichier log si n'existe pas déjà #####
if (!file.exists("log.log")) {
 log <- file("log.log", open = "w")
} else {
 log <- file("log.log", open = "a+")
}
logw <- function(textline) { writeLines(textline, log) }
#
path <- "../16_deseq2"
files <- list.files(path, pattern = "*.dge.csv", full.names = TRUE)
#
# Dans le dossier désigné, pour chaque fichier en .csv
for (file in files) {
 filename <- tools::file_path_sans_ext(basename(file))
 if (TRUE) { logw("Récupération des DGE"); 
  writeLines("Chargement des données de DGE", log) #####
  dge <- read.csv(file, header=TRUE, row.names=1)
  logw(capture.output(head(dge)))
 }
 dge <- as.data.frame(dge)
 ######################  ^ Objet DDS  v Extraction & DGE  #########################
 logw("Génération des plot conventionnels en RNAseq DGE") #####
 ########################################################################
 if (TRUE) { logw("Volcano plot");
  library(plotly)
  names(dge)[names(dge) == "log2FoldChange"] <- "log2FC"
  names(dge)[names(dge) == "padj"] <- "neglog10padj"
  dge$neglog10padj <- -log10(dge$neglog10padj)
  dge <- subset(dge, neglog10padj > 0)
  dge <- dge[!is.na(dge$neglog10padj), ]
  log2FC <- dge$log2FC
  neglog10padj <- dge$neglog10padj
  #
  volcan <- plot_ly()
  volcan <- volcan %>% layout(
                        title = "Volcano Plot",
                        xaxis = list(title = "Log2 Fold Change",
                                     titlefont = list(color = "black", family = "serif", size = 22, weight = "bold"),
                                     tickfont = list(color = "black", family = "serif", size = 22, weight = "bold")
                                    ),
                        yaxis = list(title = "-log10(Adjusted p-value)",
                                     titlefont = list(color = "black", family = "serif", size = 22, weight = "bold"),
                                     tickfont = list(color = "black", family = "serif", size = 22, weight = "bold")
                                    ),
                        showlegend = FALSE,
                        shapes = list(
                          list(type = "line", x0 = -1, y0 = 0, x1 = -1, y1 = max(dge$neglog10padj), line = list(color = "green", dash = "dash")),
                          list(type = "line", x0 = 1, y0 = 0, x1 = 1, y1 = max(dge$neglog10padj), line = list(color = "green", dash = "dash")),
                          list(type = "line", x0 = -5, y0 = 1.3, x1 = 5, y1 = 1.3, line = list(color = "green", dash = "dash"))
                        )
                       )
  threshold <- as.data.frame(subset(dge, abs(log2FC) < 1 | abs(neglog10padj) < 1.3))
  threshold_filteredlog2FC <- threshold$log2FC
  threshold_filteredneglog10padj <- threshold$neglog10padj
  threshold_label <- rownames(threshold)
  #
  underexp <- as.data.frame(subset(dge, log2FC <= -1 & neglog10padj >= 1.3))
  underexp_filteredlog2FC <- underexp$log2FC
  underexp_filteredneglog10padj <- underexp$neglog10padj
  underexp_label <- rownames(underexp)
  #
  overexp <- as.data.frame(subset(dge, log2FC >= 1 & neglog10padj >= 1.3))
  overexp_filteredlog2FC <- overexp$log2FC
  overexp_filteredneglog10padj <- overexp$neglog10padj
  overexp_label <- rownames(overexp)
  ##########
  volcan <- volcan %>% add_trace(data = threshold, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
                                 marker = list(size = 5, color = "black"),
                                 text = threshold_label
                                )
  volcan <- volcan %>% add_trace(data = underexp, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
                                 marker = list(size = 5, color = "blue"),
                                 text = underexp_label
                                )
  volcan <- volcan %>% add_trace(data = overexp, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
                                 marker = list(size = 5, color = "red"),
                                 text = overexp_label
                                )
  #
  volcan <- volcan %>% config(toImageButtonOptions = list(format = 'svg', filename = 'volcano'))
  #
  htmlwidgets::saveWidget(volcan, paste0(filename, "_", "volcano_plot.html"))
 }
 if (FALSE) { logw("Plot heatmap pheatmap"); 
  library(pheatmap)
  norm_counts <- counts(dds, normalized=TRUE)
  sig_genes <- res[which(res$padj < 0.05), ] # Filtrez les gènes significativement différentiels
  #
  library(pheatmap)
  library(tidyverse)
  library(svglite)
  library(RColorBrewer)
  # FONCTIONS
  StrSplitPointExtrNiemeElement <- function(string, n) {
   elements <- strsplit(string, "\\.")[[1]]
   if (length(elements) >= n) {
    element_n <- elements[n]
   } else {
    element_n <- NA
   }
   return(element_n)
  }
  VecNamesExtrNiemeElement <- function(vecteur, n) {
   resultats <- sapply(vecteur, StrSplitPointExtrNiemeElement, n)
   #resultats <- as.vector(resultats)
   return(resultats)
  }
  ###########
  print("Récupération des données")###
  M <- data.matrix(read.csv("Mcounts_normalized.csv", header = TRUE, row.names = 1))
  M <- M[order(rownames(M)), , drop = FALSE]
  F <- data.matrix(read.csv("Fcounts_normalized.csv", header = TRUE, row.names = 1))
  F <- F[order(rownames(F)), , drop = FALSE]
  print(head(M))
  print("Conversion en matrice de données")###
  Mrownames <- rownames(M)
  Frownames <- rownames(F)
  Crownames <- intersect(Mrownames, Frownames)
  M <- M[Crownames, , drop = FALSE]
  F <- F[Crownames, , drop = FALSE]
  Matrice <- cbind(M, F)
  #
  print(head(Matrice))
  print("Regroupement des colonnes selon Resistant ou Sensible")###
  resistants <- grep(".R.", colnames(Matrice))
  Matrice_Re <- Matrice[, resistants]
  sensibles <- grep(".S.", colnames(Matrice))
  Matrice_Se <- Matrice[, sensibles]
  Matrice <- cbind(Matrice_Re, Matrice_Se)
  #
  Matrice <- Matrice[1:100, , drop = FALSE]
  indices <- which(rownames(Matrice) %in% c(
   "HCON_00162780"
  ))
  Matrice <- Matrice[indices, ]
  #
  print(head(Matrice))
  print("Groupe(s)")###
  annotation_col <- data.frame(
   Resistance = VecNamesExtrNiemeElement(colnames(Matrice), 2),
   Sexe = VecNamesExtrNiemeElement(colnames(Matrice), 4),
   Traitement = VecNamesExtrNiemeElement(colnames(Matrice), 3)
  )
  row.names(annotation_col) <- colnames(Matrice)
  print(annotation_col)
  print("Pheatmap")###
  heatmap <- pheatmap(Matrice,
				cluster_cols = FALSE,
				cluster_rows = TRUE,
				border_col = FALSE,
				annotation_col = annotation_col
				)
  #cluster_rows = TRUE,
  ggsave(paste0(filename, "_", "heatmap.svg"), plot = heatmap, device = "svg")
  #
 }
 # Fin boucle for VVV
}
# ^^^^^^^^
EOT
#####################
# Rscript exécution #
#####################
echo "$scriptR" | $Rscript -
```

## `main.sh`

Source brute : lignes 3013–3033 de la compilation.

```bash
#!/bin/sh
#########################
# Date et heure actuelles
now() {
 echo $(date +"%Y%m%d%H%M%S%N")
}
#########################
# Création du fichier log horodaté
> log.log; chmod 777 log.log; echo "$(now)" >> log.log
############################
# Exécution des sous-scripts
for folder in $(ls -d */); do
 if [[ $folder =~ ^[0-9]*v ]]; then
  echo "$(now)" >> log.log
  echo "$folder" >> log.log
  cd $folder
  source "./sh.sh"
  cd ..
  echo "$(now)" >> log.log
 fi
done
```

## `README.md`

Source brute : lignes 3039–3083 de la compilation.

````markdown
# rnaseqdge

## Description
Pipeline d'analyse de données RNAseq paired-end : DGEA, GSEA manuel

## Installation
Préparer les dossiers, sur la même racine :
- RAW
    - fichiers de séquençage fq.gz en vrac
    - génome de référence .fa
    - annotation .gtf)
- RES
    - ici, git clone ce pipeline)
```sh
git clone https://gitlab.com/gitliro/rnaseqdge.git
```
## Usage
- Soit sur cluster Genotoul bioinfo, environnement Slurm (bientôt grâce à Singularity, sur n'importe quel cluster ayant Singularity installé)
- Soit en local
- Nommez vos fichiers selon les conditions expérimentales
par exemple R.LOC.M.1.fq.gz = échantillon Résistant, région LOC, sexe Mâle, réplicat 1

## Dépendance
- Linux
- Singularity

## Avertissements
- Pipeline non finalisé
- Images Singularity créées de FastQC à Hisat2 pour l'instant, à récupérer sur le dossier partagé dont le lien est donné, puis placer dans les dossiers étapes de pipeline correspondants
- main.sh pour exécuter les étapes successivement préparé mais non finalié, exécuter les étapes (.sh), par dossier, manuellement
- En raison d'une panne du cluster HPC INRAe dernièrement, des scirpts ont été modifiés manuellement au lieu de commit. Ainsi, en cas de bug à cause des retours chariot \r Windows, passer les scripts problématiques dans
```sh
sed -i 's/\r//' LeScriptProblematique.sh
```

## Licence
Sous licence CeCILL. Appartient à INRAe UMR 1436.
Citer l'auteur. Mentionner licence CeCILL. Si modifications garder la traçabilité en citant le code source. Contagion de la licence si redistribution. Copyleft, ne pas intégrer à logiciel propriétaire.

## Auteurs
- LIOUTAUD Robin
- INRAe Centre Occitanie-Toulouse, UMR INRAe/ENVT 1436 InTheRes, Toulouse, France

## Contact
UMR 1436 INRAe
````
