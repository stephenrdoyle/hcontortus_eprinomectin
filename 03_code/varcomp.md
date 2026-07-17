# Documentation of the comparative variant analysis pipeline

> **Purpose of this document.** This page converts a raw laboratory notebook into scientific documentation intended to accompany a Git repository associated with a publication. It describes the **commands, parameters, tools, files, filters, calculations and visualisations actually recorded in the supplied code**.
>
> The purpose is not to present the repository as a fully reproducible or immediately reanalysable workflow. Rather, it is to allow readers to determine precisely **which operations were programmed and according to which rules**.

## Summary

This repository implements a comparative variant analysis between two groups of samples:

- **resistant** samples, identified by a name beginning with `R_` or, after renaming, by a column beginning with `R`;
- **susceptible** samples, identified by a name beginning with `S_` or by a column beginning with `S`.

The pipeline assumes that the variants have already been:

1. called individually, notably with `bcftools mpileup/call`, according to the script help text;
2. annotated with **SnpEff**;
3. compressed as `VCF.gz` files;
4. indexed with a `.tbi` file;
5. generated against the same reference genome.

The pipeline then performs the following operations:

```text
Individual SnpEff-annotated VCF.gz files
        │
        ├── relabelling of sample names
        ├── multisample merging with bcftools
        ├── conversion to TSV
        ├── calculation of AF from AD and DP
        ├── differential R/S filtering
        ├── assignment to genes and exons from the GTF
        ├── extraction of SnpEff effects
        └── visualisations:
                ├── sliding windows
                ├── individual variants
                ├── interactive gene-level aggregation
                └── static gene-level aggregation
```

## Exact scope

The repository **does not contain** the complete commands used to generate the initial individual VCF files. Variant calling and SnpEff annotation are external prerequisites. The supplied code effectively starts from individual VCF files that have already been annotated and indexed.

Two different methods are used to link variants, genes and effects:

| Script | Function |
|---|---|
| `4_gtf_vcf_join.r` | spatial overlap between variant positions and gene/exon intervals from a simplified GTF |
| `4_variants_to_gene_effect_table.r` | extraction of SnpEff `ANN` annotations, association with `HCON_*` identifiers, and subsequent joining to genes from the GTF |

The two scripts generate different tables and should be regarded as either complementary outputs or competing approaches, rather than as a single, perfectly equivalent redundant step.

# Directory structure

```text
.
├── 0_raw_and_do_not_use/
├── 0_refseq/
│   └── <reference.gtf>
├── 0_Rscripts/
│   ├── 4_gtf_vcf_join.r
│   ├── 4_variants_to_gene_effect_table.r
│   ├── 5_var_density.r
│   ├── 6_plotly_var.r
│   ├── 7_ggplot_gene.r
│   └── 7_plotly_gene.r
├── 0_samples/
│   ├── R_<sample>/
│   │   ├── <sample>.vcf.gz
│   │   └── <sample>.vcf.gz.tbi
│   └── S_<sample>/
│       ├── <sample>.vcf.gz
│       └── <sample>.vcf.gz.tbi
├── 0_tools/
├── README.md
└── VarEffPlot.sh
```

Directories numbered from `1_` to `7_` are created automatically during execution.

# Inputs and prerequisites

## Individual VCF files

Each subdirectory within `0_samples/` must contain:

```text
<sample>.vcf.gz
<sample>.vcf.gz.tbi
```

The script selects the **first file** matching `*.vcf.gz` in each directory:

```bash
vcf="$(find "$sample" -maxdepth 1 -type f -name "*.vcf.gz" | head -n 1)"
```

Each directory name must begin with:

| Prefix | Group |
|---|---|
| `R_` | resistant |
| `S_` | susceptible |

The directory name becomes the new sample name in the VCF.

## GTF reference

The first `*.gtf` file found directly within `0_refseq/` is used:

```bash
find ./0_refseq -maxdepth 1 -type f -name '*.gtf' | head -n 1
```

The pipeline assumes, in particular, that:

- chromosomes in the GTF and VCF use the same identifiers;
- gene attributes follow the syntax `gene_id "..."`;
- the expected identifiers may follow the pattern `HCON_\d+`.

## Genotype-field format

The allele-frequency calculation implicitly assumes that the sample fields are arranged as follows after splitting on `:`:

```text
field 1: genotype or another field
field 2: another field
field 3: DP
field 4: AD
```

The code uses:

```awk
dp = a[3]
ad = a[4]
```

However, the exact order specified in the `FORMAT` field is not read dynamically. The pipeline therefore depends on a specific VCF format.

# Software and dependencies

No software version is explicitly pinned in this compilation.

| Tool or package | Use |
|---|---|
| Bash | overall orchestration |
| bcftools | reading sample names, relabelling, indexing and merging VCF files |
| SnpEff | prior annotation of consequences in the `ANN` field |
| awk | conversion, AF calculation, filtering and initial association with the GTF |
| gunzip | decompression of the merged VCF |
| curl or wget | download of auxiliary `.tar` files |
| R / Rscript | tabular processing and figure generation |
| data.table | rapid input, interval joins and aggregation |
| stringr | extraction of identifiers and annotations |
| readr | reading and writing TSV files |
| dplyr | transformations and aggregations |
| ggplot2 | static plots |
| ggrepel | gene labelling |
| viridis | density colour scale |
| plotly | interactive plots |
| htmlwidgets | HTML export |
| rmarkdown | Pandoc verification and standalone HTML export |
| Pandoc | generation of `selfcontained` HTML files |

The script also downloads:

```text
bcftools.tar
r-plotly.tar
r-ggplot2.tar
```

from an IONOS URL. However, these archives are **neither extracted nor added to the `PATH` or the R library path** in the supplied script. In practice, execution therefore relies on installations that are already accessible in the environment, unless an undocumented external intervention is performed.

# Overall orchestration

## Help

```bash
bash VarEffPlot.sh -h
```

Options available for the static gene-level figure:

| Option | Effect |
|---|---|
| `-l <bp>` | left genomic boundary |
| `-r <bp>` | right genomic boundary |
| `-b <value>` | lower bound for `AF_R - AF_S` |
| `-n` | add gene identifiers |
| `-f` | force the Y-axis to range from 0 to 1 |

## Example command

```bash
bash VarEffPlot.sh -l 1000000 -r 5000000 -b 0.2 -n
```

## Restart logic

Each step is guarded by a test for the existence of its output directory:

```bash
if [ ! -d "4_join" ]; then
    ...
fi
```

Consequently:

- if the directory exists, the step is skipped;
- to recompute a step, the README recommends deleting the corresponding directory and rerunning the pipeline;
- the mere existence of the directory, even if it is empty or incomplete, is sufficient to prevent re-execution.

# Detailed steps

## 0. Downloading auxiliary files

```bash
mkdir -p "0_tools"
```

Function:

```bash
download_tool() {
    file=$1
    url="https://storage.ionos.fr/s/cQNq7E5iQ763Rcw/download?path=%2F&files=$file"
    out="./0_tools/$file"

    if [ ! -s "$out" ]; then
        if command -v curl >/dev/null 2>&1; then
            curl -fL "$url" -o "$out"
        elif command -v wget >/dev/null 2>&1; then
            wget -O "$out" "$url"
        else
            exit 1
        fi
    fi
}
```

Requested files:

```bash
download_tool bcftools.tar
download_tool r-plotly.tar
download_tool r-ggplot2.tar
```

## 1. Relabelling individual VCF files

**Output directory:** `1_vcf_reheadered/`

For each sample directory:

```bash
samplename=$(basename "$sample")
oldname="$(bcftools query -l "$vcf")"

printf "%s\t%s\n" "$oldname" "$samplename" > "$mapfile_tmp"

bcftools reheader \
    -s "$mapfile_tmp" \
    -o "./1_vcf_reheadered/$samplename/$samplename.vcf.gz" \
    "$vcf"

bcftools index -t \
    "./1_vcf_reheadered/$samplename/$samplename.vcf.gz"
```

The script implicitly assumes that each VCF contains a single sample name. If `bcftools query -l` returns several lines, the renaming table constructed by `printf` will be incorrect.

## 2. Merging VCF files

**Output directory:** `2_vcf_merged/`

```bash
bcftools merge \
    -Oz \
    -o ./2_vcf_merged/all.vcf.gz \
    $(find ./1_vcf_reheadered/ -type f -name '*.vcf.gz' -print0 |
      xargs -0 printf '%s ')

bcftools index -t ./2_vcf_merged/all.vcf.gz
```

The merge is performed without additional options defining the handling of multiallelic alleles, absent sites or missing genotypes. The defaults of the installed `bcftools` version therefore apply.

## 3A. GTF simplification

**Output directory:** `3_gtf_filtered/`  
**File:** `gtf.tsv`

Only `gene` and `exon` records are retained.

```bash
awk -F'\t' 'BEGIN{OFS="\t"}
  $0 ~ /^#/ {next}
  ($3!="gene" && $3!="exon") {next}
  {
    chrom=$1
    start=$4
    end=$5

    gid="NA"
    if(match($0, /gene_id "[^"]+"/)){
      gid=substr($0, RSTART+9, RLENGTH-10)
    }

    t=($3=="gene" ? "g" : "e")
    print chrom, gid, t, start, end
  }' "$GTF" > "./3_gtf_filtered/gtf.tsv"
```

Output format:

```text
chromosome    gene_id    type    start    end
```

where:

- `g` denotes a gene record;
- `e` denotes an exon record.

## 3B. Conversion and filtering of the merged VCF

**Output directory:** `3_vcf_filtered/`  
**Final file:** `vcf.tsv`

### 3B.1. Removal of `##` metadata lines

```bash
gunzip -c "$cur" |
    awk 'BEGIN{OFS="\t"} /^##/ {next} {print}' \
    > 01.no_double_hash.tsv
```

The `#CHROM` line is retained.

### 3B.2. Column selection

The script retains:

```text
#CHROM
POS
REF
ALT
INFO
sample columns
```

The following columns are removed:

```text
ID
QUAL
FILTER
FORMAT
```

Removal of `FORMAT` is important: the sample fields are subsequently interpreted according to a fixed order, without retaining the description of that order.

### 3B.3. Counting SnpEff impact categories

The number of occurrences is calculated directly within `INFO`:

```awk
high     = gsub(/\|HIGH\|/, "", tmp)
moderate = gsub(/\|MODERATE\|/, "", tmp)
modifier = gsub(/\|MODIFIER\|/, "", tmp)
low      = gsub(/\|LOW\|/, "", tmp)
```

Added columns:

```text
HIGH
MODERATE
MODIFIER
LOW
```

These values represent the number of corresponding substrings in `INFO`, not necessarily the number of biologically independent effects.

### 3B.4. Calculation of allele frequencies

For each sample:

```awk
dp = a[3]
ad = a[4]
split(ad, adp, ",")
altad = (length(adp)>=2 ? adp[2] : adp[1])
AF = altad / dp
```

The value is formatted to two decimal places:

```awk
sprintf("%.2f", AF)
```

Missing or invalid cases:

```text
./.  → 0.00
.    → 0.00
DP ≤ 0 → 0.00
missing or invalid AD → 0.00
```

Missing data are therefore converted to an allele frequency of zero, which is not biologically equivalent to a confident observation of the reference allele.

### 3B.5. Depth

The minimum and maximum `DP` values across samples are calculated:

```text
minDP
maxDP
```

The following columns are created but systematically remain `NA`:

```text
minQUAL
maxQUAL
```

Because the VCF `QUAL` column was removed earlier, no variant-quality metric is actually summarised.

### 3B.6. Filtering susceptible samples

A variant is removed if its allele frequency exceeds `0.1` in **at least one** sample whose name begins with `S`:

```awk
if($c!="NA" && ($c+0)>0.1){ bad=1 }
```

The threshold is strict:

```text
AF > 0.1
```

A value equal to `0.1` is retained.

### 3B.7. Filtering resistant samples

Parameter:

```bash
MIN_R_COUNT=1
```

A variant is retained if:

```text
AF > 0.4
```

in at least `MIN_R_COUNT` columns beginning with `R`.

In the supplied configuration:

```text
at least one resistant sample with AF > 0.4
```

A value equal to `0.4` does not pass the filter.

### 3B.8. Preliminary gene assignment

An `awk` step adds an `HCON_ID` column by identifying genes from the GTF that contain the variant position.

```text
CHROM    HCON_ID    POS    REF    ALT    INFO    ...
```

Multiple genes are concatenated with commas. Variants outside genes are assigned `NA`.

This column is added to the final `3_vcf_filtered/vcf.tsv` file, but the spatial-overlap R script subsequently performs its own variant–gene assignment from the simplified GTF.

## 4A. Spatial VCF–GTF join

**Script:** `0_Rscripts/4_gtf_vcf_join.r`  
**Output:** `4_join/var.tsv`

Packages:

```r
library(data.table)
library(stringr)
```

### Input and checks

Required VCF columns:

```r
c("CHROM", "POS", "REF", "ALT", "INFO")
```

The simplified GTF is assigned the following column names:

```r
c("chromosome", "gene_id", "type", "start", "end")
```

### Gene overlap

Each variant is converted into a point interval:

```r
vcf_pos <- vcf[, .(
    row_id__,
    chromosome = CHROM,
    vstart = POS,
    vend = POS
)]
```

Join:

```r
jg <- foverlaps(
    vcf_pos,
    genes_iv,
    by.x = c("chromosome", "vstart", "vend"),
    by.y = c("chromosome", "start", "end"),
    type = "within",
    nomatch = NA
)
```

### Exon detection

A second join tests whether the position lies within an exon of the same gene.

Region code:

| Value | Interpretation |
|---|---|
| `e` | variant within an exon |
| `g` | variant within a gene, but with no detected exon overlap |
| empty string | variant not assigned to a gene |

The code therefore does not explicitly distinguish intronic, UTR, promoter or intergenic locations in `region_code`.

### Gene identifiers in `MODIFIER` annotations

The `extract_modifier_gene_ids()` function:

1. extracts `ANN=` from `INFO`;
2. splits annotations on commas;
3. retains only annotations whose third field is `MODIFIER`;
4. searches fields 4 and 5 for patterns matching `HCON_\d+`;
5. concatenates unique identifiers using `$`.

Output:

```text
modifier_gene_ids
```

### Leading columns

```text
chromosome
gene
gene_start
gene_end
region_code
variant_pos
ref
alt
modifier_gene_ids
```

The remaining columns from the VCF-derived TSV are then retained.

## 4B. Variant–gene–effect table

**Script:** `0_Rscripts/4_variants_to_gene_effect_table.r`  
**Main output:** `4_gene_var_eff/outnew`

Intermediate outputs:

```text
4_gene_var_eff/prepared.gtf.csv
4_gene_var_eff/prepared.vcf.csv
```

### GTF preparation

Only records satisfying `feature == "gene"` are retained.

Extracted attributes:

```text
gene_id
gene_name
```

The `gene:` prefix is removed from `gene_id`. If `gene_name` is absent, `gene_id` is used.

### Extraction of the SnpEff `ANN` field

The `ANN` field is split according to the 16 expected standard fields:

```text
allele
annotation
impact
gene_name
gene_id
feature_type
feature_id
transcript_biotype
rank
hgvs_c
hgvs_p
cdna_pos
cds_pos
aa_pos
distance
warnings
```

The script extracts all patterns matching:

```r
HCON_\d+
```

from `gene_id` and `gene_name`.

One row is created for each combination:

```text
variant × SnpEff effect × HCON_ID
```

Identifiers absent from the GTF are removed.

### Final table

Columns:

```text
chromosome
gene
gene_start
gene_end
variant_pos
variant_ref
variant_alt
<sample columns>
snpeff_effect
snpeff_impact
```

This output does not appear to be used by the subsequent plotting scripts, which primarily read `4_join/var.tsv`.

## 5. Sliding windows and variant density

**Script:** `5_var_density.r`  
**Input:** `4_join/var.tsv`  
**Output:** `5_var_density/sliding_AF_R_minus_S_by_chr_plotly.html`

Packages:

```r
data.table
ggplot2
viridis
plotly
htmlwidgets
```

Effective parameters:

```r
window_size <- 500000
step_size   <- 10000
```

The header comment states a step size of `100000 bp`, whereas the executed variable uses `10000 bp`.

### Group detection

```r
s_cols <- grep("^S_", names(dt), value = TRUE)
r_cols <- grep("^R_", names(dt), value = TRUE)
```

Unlike some of the other scripts, which use `^S` and `^R`, this script specifically requires names beginning with `S_` and `R_`.

### Per-variant calculation

```r
mean_AF_S = mean of susceptible-sample columns
mean_AF_R = mean of resistant-sample columns
delta_AF  = mean_AF_R - mean_AF_S
```

### Window-level aggregation

For each chromosome:

```r
starts <- seq(1, max_pos, by = step_size)
win_end <- win_start + window_size - 1
```

The windows overlap because:

```text
window = 500 kb
step = 10 kb
```

Metrics:

```text
n_variants
mean_delta
```

Empty windows are assigned `n_variants = 0` and `mean_delta = NA`.

### Figure

- one panel per chromosome;
- black `mean_delta` curve;
- horizontal line at zero;
- trapezoids beneath the curve coloured by `n_variants`;
- standalone Plotly export.

## 6. Interactive variant-level visualisation

**Script:** `6_plotly_var.r`  
**Input:** `4_join/var.tsv`

Outputs:

```text
6_plotly_var/<chromosome>.html
6_plotly_var/top10pct_variants_by_absDiffAF.tsv
```

### AF columns

```r
af_cols   <- grep("^[RS]", names(df), value = TRUE)
r_samples <- grep("^R", names(df), value = TRUE)
s_samples <- grep("^S", names(df), value = TRUE)
```

This rule may select any column beginning with `R` or `S`, rather than only sample columns, if other column names share these initial letters.

### Primary impact

If the counters are present, the priority order is:

```text
HIGH > MODERATE > MODIFIER > LOW > NONE
```

### Frequency difference

```r
mean_R  <- mean resistant-sample AF
mean_S  <- mean susceptible-sample AF
diff_AF <- mean_R - mean_S
```

### Top 10%

Threshold:

```r
quantile(abs(diff_AF), probs = 0.90)
```

Variants with `|diff_AF|` greater than or equal to the 90th percentile are exported.

### Figure

Each bubble corresponds to one row of the table and may therefore represent a variant–gene combination when the join has generated multiple rows for a single variant.

Axes:

```text
X = variant position
Y = mean(AF_R) - mean(AF_S)
```

Colour:

```text
primary impact
```

Background:

- blue rectangles for gene intervals;
- green rectangles for exons when their boundaries are available;
- otherwise, point marks at positions annotated as `e`.

A standalone HTML file is generated for each chromosome. Pandoc is required.

## 7A. Static gene-level aggregated visualisation

**Script:** `7_ggplot_gene.r`  
**Input:** `4_join/var.tsv`

Outputs:

```text
7_ggplot_gene/ggplot_<chromosome>.pdf
7_ggplot_gene/top10pct_genes_by_absDiffAF.tsv
```

Only variants satisfying:

```r
region_code %in% c("g", "e")
```

are retained.

### Aggregation

Grouping variables:

```text
chromosome
gene
gene_start
gene_end
```

Metrics:

```text
gene_mid
gene_len_bp
n_variants
mean AF per sample
best_impact
mean_R
mean_S
diff_AF
variants_per_kb
```

Density:

```r
variants_per_kb = n_variants * 1000 / gene_len_bp
```

The number of unique variants is based on:

```text
position:REF:ALT
```

### Top 10% of genes

```r
quantile(abs(gene_summary$diff_AF), probs = 0.90)
```

### Figure

```text
X = gene midpoint in Mb
Y = mean AF_R - mean AF_S
size = variants per kb
colour = highest-priority impact
```

In the supplied version, however, all impact categories use the same hexadecimal colour:

```r
"#440154"
```

The comments to the right of the colour definitions suggest that other colours had been considered, but they are not active.

Gene labels are added with `ggrepel` when the `-n` option is supplied.

## 7B. Interactive gene-level aggregated visualisation

**Script:** `7_plotly_gene.r`  
**Input:** `4_join/var.tsv`

Expected outputs are written to:

```text
7_plotly_gene/
```

This script reproduces the gene-level aggregation but adds more detailed extraction from `INFO/ANN`.

### Selection of a SnpEff annotation

For each variant–gene combination, `extract_best_annotation_fast()` evaluates annotations according to:

1. gene match;
2. alternative-allele match;
3. `protein_coding` biotype;
4. presence of a protein HGVS annotation;
5. effect priority;
6. SnpEff impact level.

Representative effect scores include:

| Effect | Score |
|---|---:|
| `stop_gained` | 100 |
| `frameshift_variant` | 95 |
| `start_lost` | 90 |
| canonical splice site | 88 |
| in-frame insertion/deletion | 80 |
| `missense_variant` | 70 |
| synonymous | 60 |
| splice region | 50 |
| UTR | 39–40 |
| intron | 30 |
| non-coding | 20 |
| upstream/downstream | 9–10 |
| intergenic | 5 |

If a protein HGVS annotation is available, the `p.` prefix is removed and the resulting notation becomes the displayed effect label.

### Gene-level aggregation

The script calculates, among other quantities:

```text
number of variants
gene length
variants per kb
mean AF per sample
mean R AF
mean S AF
diff_AF
primary impact
list of retained variants and effects
```

The interactive figures allow information to be inspected by hovering.

# Main outputs

| Directory | File or output type |
|---|---|
| `1_vcf_reheadered/` | renamed individual VCF files and `.tbi` indices |
| `2_vcf_merged/` | `all.vcf.gz`, `all.vcf.gz.tbi` |
| `3_gtf_filtered/` | `gtf.tsv` |
| `3_vcf_filtered/` | intermediate steps `01` to `06`, followed by `vcf.tsv` |
| `4_join/` | `var.tsv` |
| `4_gene_var_eff/` | `prepared.gtf.csv`, `prepared.vcf.csv`, `outnew` |
| `5_var_density/` | sliding-window HTML figure |
| `6_plotly_var/` | per-chromosome HTML files and top 10% of variants |
| `7_plotly_gene/` | per-chromosome HTML figures |
| `7_ggplot_gene/` | per-chromosome PDF files and top 10% of genes |

# Critical thresholds and parameters

| Parameter | Value |
|---|---:|
| maximum permitted susceptible-sample AF | no S sample with `AF > 0.1` |
| minimum resistant-sample AF | at least one R sample with `AF > 0.4` |
| `MIN_R_COUNT` | 1 |
| exported AF precision | 2 decimal places |
| sliding-window size | 500,000 bp |
| effective step size | 10,000 bp |
| top variants | upper decile of `|AF_R - AF_S|` |
| top genes | upper decile of `|AF_R - AF_S|` |

# Points requiring caution

1. **Software versions are absent.** No versions are recorded for bcftools, SnpEff, R or the R packages.
2. **Variant calling is external.** The `mpileup/call` commands, their parameters and any preceding filters are not included in the repository.
3. **SnpEff annotation is external.** The SnpEff database and its version are not specified.
4. **Archives are downloaded but not installed.** The `.tar` archives for bcftools and the R packages are neither extracted nor explicitly used.
5. **Only the first file is selected.** The first GTF and the first VCF in each directory are selected without checking for ambiguity.
6. **Single-sample VCF assumption.** Relabelling does not explicitly handle multiple sample names in a source file.
7. **`FORMAT` is removed.** Calculation of DP and AD relies on fields 3 and 4, rather than on keys read from the `FORMAT` field.
8. **Missing data are converted to AF = 0.** This convention may artificially increase the R/S contrast.
9. **AF is rounded before filtering and analysis.** Precision is reduced to two decimal places.
10. **`QUAL` is not used.** `minQUAL` and `maxQUAL` remain `NA`.
11. **No explicit DP filter is applied.** `minDP` and `maxDP` are calculated, but no depth threshold is enforced.
12. **Asymmetric filters.** A single R sample above 0.4 is sufficient, whereas a single S sample above 0.1 excludes the variant.
13. **Inconsistent prefixes.** Some scripts search for `R_`/`S_`, whereas others search only for `R`/`S`.
14. **A Cartesian join is possible.** A variant overlapping several genes produces several rows.
15. **Broad `g` code.** Any position within a gene but outside a detected exon is assigned `g`, without distinguishing intronic from UTR sequence.
16. **Two association methods.** `4_join/var.tsv` and `4_gene_var_eff/outnew` are not strictly equivalent.
17. **Impact categories are counted by text pattern.** `HIGH`, `MODERATE`, and related values are counts of text occurrences within `INFO`.
18. **Windows overlap extensively.** A 500 kb window with a 10 kb step induces strong visual autocorrelation.
19. **Comment–code inconsistency.** The comment states a 100 kb step, whereas the code uses 10 kb.
20. **The top 10% is based on available rows.** Ties at the quantile threshold may result in more than 10% of rows being selected.
21. **Variant-level figures may contain duplicates.** A multi-gene join may display several bubbles for the same variant.
22. **Identical colours in the static plot.** In the supplied version, the gene-level ggplot does not visually distinguish impact categories.
23. **Restart behaviour depends only on directories.** An incomplete directory prevents automatic re-execution.
24. **No global logging.** Versions, dates, standard output/error and exit codes are not recorded centrally.
25. **No chromosome-compatibility check.** Differences in chromosome naming between the GTF and VCF result in unassigned variants.

# Information to add for the publication repository

| Item | Status |
|---|---|
| bcftools version | to be provided |
| SnpEff version and database used | to be provided |
| R and package versions | to be provided if historically available |
| reference genome and annotation | accession and version to be provided |
| upstream variant-calling commands | to be documented separately |
| filters applied before this pipeline | to be documented |
| exact biological definition of the R and S groups | to be added |
| sample-to-condition table | to be added |
| depth and quality criteria | to be specified, or their absence stated explicitly |
| commit associated with the published figures | to be fixed |
| execution date | to be added |
| data DOI and accession | to be added |

# Appendix — files reproduced without modification

The files below are reproduced exactly as they appear in the compilation. `.gitkeep` files containing only the indication `[Binary content omitted]` are not reproduced as analytical code.


## `0_Rscripts/4_gtf_vcf_join.r`

Raw source: lines 19–145 of the compilation.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

vcf_file <- "./3_vcf_filtered/vcf.tsv"
gtf_file <- "./3_gtf_filtered/gtf.tsv"
out_dir  <- "./4_join"
out_file <- file.path(out_dir, "var.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

extract_modifier_gene_ids <- function(info_vec) {
  out <- rep("", length(info_vec))
  for (i in seq_along(info_vec)) {
    info <- info_vec[i]
    if (is.na(info) || !grepl("(^|;)ANN=", info)) next
    ann_raw <- sub(".*(?:^|;)ANN=", "", info)
    ann_raw <- sub(";.*", "", ann_raw)
    if (ann_raw == "") next

    anns <- strsplit(ann_raw, ",", fixed = TRUE)[[1]]
    genes <- character(0)

    for (ann in anns) {
      f <- strsplit(ann, "\\|", perl = TRUE)[[1]]
      if (length(f) < 5) next
      if (f[3] != "MODIFIER") next
      ids <- str_extract_all(paste(f[4], f[5], sep=";"), "HCON_\\d+")[[1]]
      if (length(ids)) genes <- c(genes, ids)
    }
    genes <- unique(genes)
    out[i] <- if (length(genes)) paste(genes, collapse="$") else ""
  }
  out
}

# --- VCF TSV ---
vcf <- fread(vcf_file, sep="\t", header=TRUE, na.strings=c("", ".", "NA"))
if ("#CHROM" %in% names(vcf)) setnames(vcf, "#CHROM", "CHROM")

need <- c("CHROM","POS","REF","ALT","INFO")
miss <- setdiff(need, names(vcf))
if (length(miss)) stop("Colonnes manquantes dans vcf.tsv: ", paste(miss, collapse=", "))

vcf[, POS := as.integer(POS)]
vcf[, row_id__ := .I]

# --- GTF TSV simplifié (5 colonnes) ---
gtf <- fread(gtf_file, sep="\t", header=FALSE,
             col.names=c("chromosome","gene_id","type","start","end"))

gtf[, start := as.integer(start)]
gtf[, end   := as.integer(end)]
gtf <- gtf[!is.na(start) & !is.na(end)]

gtf_genes <- gtf[type == "g" & gene_id != "" & gene_id != "NA",
                 .(chromosome, gene_id, gene_start=start, gene_end=end)]
gtf_exons <- gtf[type == "e" & gene_id != "" & gene_id != "NA",
                 .(chromosome, gene_id, exon_start=start, exon_end=end)]

if (nrow(gtf_genes) == 0) stop("Aucun gène (type=g) dans ", gtf_file)

# --- Intervalles variants ---
vcf_pos <- vcf[, .(row_id__, chromosome=CHROM, vstart=POS, vend=POS)]
setkey(vcf_pos, chromosome, vstart, vend)

# --- Join gènes ---
genes_iv <- copy(gtf_genes)
setnames(genes_iv, c("gene_start","gene_end"), c("start","end"))
setkey(genes_iv, chromosome, start, end)

jg <- foverlaps(vcf_pos, genes_iv,
                by.x=c("chromosome","vstart","vend"),
                by.y=c("chromosome","start","end"),
                type="within", nomatch=NA)

jg <- jg[, .(row_id__, gene_id, gene_start=start, gene_end=end)]

# --- Join exons ---
if (nrow(gtf_exons) > 0) {
  exons_iv <- copy(gtf_exons)
  setnames(exons_iv, c("exon_start","exon_end"), c("start","end"))
  setkey(exons_iv, chromosome, start, end)

  je <- foverlaps(vcf_pos, exons_iv,
                  by.x=c("chromosome","vstart","vend"),
                  by.y=c("chromosome","start","end"),
                  type="within", nomatch=NA)
  exon_hits <- unique(je[!is.na(gene_id), .(row_id__, gene_id)])
  exon_hits[, in_exon__ := TRUE]
} else {
  exon_hits <- data.table(row_id__=integer(), gene_id=character(), in_exon__=logical())
}

# --- Assemblage ---
res <- merge(jg, vcf, by="row_id__", all.y=TRUE, allow.cartesian=TRUE)
res <- merge(res, exon_hits, by=c("row_id__","gene_id"), all.x=TRUE)

res[, region_code := ""]
res[!is.na(gene_id), region_code := "g"]
res[!is.na(in_exon__) & in_exon__ == TRUE, region_code := "e"]

mod_tbl <- vcf[, .(row_id__, modifier_gene_ids = extract_modifier_gene_ids(INFO))]
res <- merge(res, mod_tbl, by="row_id__", all.x=TRUE)

setnames(res, c("CHROM","POS","REF","ALT"),
             c("chromosome","variant_pos","ref","alt"),
             skip_absent=TRUE)

res[, gene := fifelse(is.na(gene_id), "", gene_id)]
res[is.na(modifier_gene_ids), modifier_gene_ids := ""]
res[is.na(region_code), region_code := ""]

front <- c("chromosome","gene","gene_start","gene_end","region_code",
           "variant_pos","ref","alt","modifier_gene_ids")

vcf_other <- setdiff(names(vcf), c("row_id__","CHROM","POS","REF","ALT"))
keep <- c(front, intersect(vcf_other, names(res)))
keep <- setdiff(unique(keep), c("row_id__","gene_id","in_exon__"))

final <- res[, ..keep]
setorderv(final, c("chromosome","variant_pos","gene"), na.last=TRUE)

fwrite(final, out_file, sep="\t", na="")
cat("OK -> ", out_file, "\n", sep="")
```

## `0_Rscripts/4_variants_to_gene_effect_table.r`

Raw source: lines 151–313 of the compilation.

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

var_file  <- "./3_vcf_filtered/vcf.tsv"
gtf_file <- system("find ./0_refseq -maxdepth 1 -type f -name '*.gtf' | head -n 1", intern = TRUE)
out_final <- "./4_gene_var_eff/outnew"

gtf_out_file <- "./4_gene_var_eff/prepared.gtf.csv"
vcf_out_file <- "./4_gene_var_eff/prepared.vcf.csv"

## ------------------------------------------------------------------------
## 1) PREP GTF  -> prepared.gtf.csv
## ------------------------------------------------------------------------

gtf_lines <- readLines(gtf_file)
gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]

parse_gtf_line <- function(line) {
  if (grepl("\t", line, fixed = TRUE)) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  } else {
    parts <- strsplit(line, "\\s+")[[1]]
  }
  if (length(parts) < 9) parts <- c(parts, rep("", 9 - length(parts)))
  fixed <- parts[1:8]
  if (length(parts) > 9) {
    attribute <- paste(parts[9:length(parts)], collapse = " ")
  } else {
    attribute <- parts[9]
  }
  c(fixed, attribute)
}

gtf_mat   <- t(vapply(gtf_lines, parse_gtf_line, FUN.VALUE = character(9)))
gtf      <- as.data.table(gtf_mat)
setnames(gtf, c("seqname","source","feature","start","end",
                "score","strand","frame","attribute"))

gtf[, start := suppressWarnings(as.integer(start))]
gtf[, end   := suppressWarnings(as.integer(end))]

gtf_genes <- gtf[feature == "gene"]

extract_attr <- function(x, key) {
  m <- stringr::str_match(x, paste0(key, " \"([^\"]+)\""))
  m[, 2]
}

gtf_genes[, gene_id   := extract_attr(attribute, "gene_id")]
gtf_genes[, gene_name := extract_attr(attribute, "gene_name")]
gtf_genes[, gene_id   := sub("^gene:", "", gene_id)]
gtf_genes[is.na(gene_name) | gene_name == "", gene_name := gene_id]

gtf_prep <- gtf_genes[, .(
  gene_id,
  gene_name,
  gene_start = start,
  gene_end   = end,
  chromosome_gtf = seqname
)]

fwrite(gtf_prep, gtf_out_file, sep = ",", na = "")
valid_gene_ids <- unique(gtf_prep$gene_id)

## ------------------------------------------------------------------------
## 2) PREP VARIANTS/ANN -> prepared.vcf.csv
## ------------------------------------------------------------------------

vars <- fread(var_file)

if ("#CHROM" %in% names(vars)) setnames(vars, "#CHROM", "CHROM")

req <- c("CHROM","POS","REF","ALT","INFO")
miss <- setdiff(req, names(vars))
if (length(miss) > 0) {
  stop(paste("Colonnes manquantes dans le fichier de variants:", paste(miss, collapse = ", ")))
}

## repérer les colonnes "samples" (RARA, RBET, RBUN, RMOU, SCHI, SLUC, etc.)
## = toutes les colonnes qui ne sont pas CHROM, POS, REF, ALT, INFO, ANN_raw
sample_cols <- setdiff(names(vars), c("CHROM","POS","REF","ALT","INFO","ANN_raw"))

vars[, ANN_raw := sub(".*ANN=", "", INFO)]
vars[, ANN_raw := sub(";.*", "", ANN_raw)]
vars[ANN_raw == INFO | ANN_raw == "", ANN_raw := NA_character_]

# on ne garde que ceux avec ANN
ann_dt <- vars[!is.na(ANN_raw), c("CHROM","POS","REF","ALT", sample_cols, "ANN_raw"), with = FALSE]

# éclater ANN sur les virgules -> 1 ligne / effet
ann_dt[, ANN_list := strsplit(ANN_raw, ",", fixed = TRUE)]
ann_exp <- ann_dt[, .(ANN = unlist(ANN_list)),
                  by = c("CHROM","POS","REF","ALT", sample_cols)]

# séparer le champ ANN selon le format SnpEff
fields <- c("allele","annotation","impact","gene_name","gene_id",
            "feature_type","feature_id","transcript_biotype",
            "rank","hgvs_c","hgvs_p",
            "cdna_pos","cds_pos","aa_pos","distance","warnings")

ann_exp[, c(fields) := tstrsplit(ANN, "\\|", fixed = FALSE, fill = "")]

# extraire TOUS les HCON_ids possibles (gene_id + gene_name)
ann_exp[, source_ids   := paste(gene_id, gene_name, sep = ";")]
ann_exp[, ids_list     := stringr::str_extract_all(source_ids, "HCON_\\d+")]

# garder seulement les lignes avec au moins un HCON et dupliquer
ann_exp2 <- ann_exp[lengths(ids_list) > 0]

## construire vcf_prep : une ligne par combinaison
## (variant, effet, HCON_id), en gardant aussi les colonnes samples
by_cols <- c("CHROM","POS","REF","ALT", sample_cols, "annotation","impact")

vcf_prep <- ann_exp2[
  ,
  .(gene_id = unlist(ids_list)),
  by = by_cols
]

# ne garder que les gènes réellement présents dans le GTF
vcf_prep <- vcf_prep[gene_id %in% valid_gene_ids]

# renommer les colonnes de base
setnames(vcf_prep,
         old = c("CHROM","POS","REF","ALT","annotation","impact"),
         new = c("chromosome","variant_pos","variant_ref","variant_alt",
                 "snpeff_effect","snpeff_impact"))

fwrite(vcf_prep, vcf_out_file, sep = ",", na = "")

## ------------------------------------------------------------------------
## 3) JOINTURE prepared.vcf.csv + prepared.gtf.csv -> table finale
## ------------------------------------------------------------------------

setkey(gtf_prep, gene_id)
setkey(vcf_prep, gene_id)

final <- merge(vcf_prep, gtf_prep, by = "gene_id", all.x = TRUE)

## ordre des colonnes dans la sortie finale
result_cols <- c(
  "chromosome",
  "gene_id",
  "gene_start",
  "gene_end",
  "variant_pos",
  "variant_ref",
  "variant_alt",
  sample_cols,
  "snpeff_effect",
  "snpeff_impact"
)

result <- final[, ..result_cols]
setnames(result, "gene_id", "gene")

setorder(result, chromosome, gene, variant_pos, snpeff_effect)

fwrite(result, out_final, sep = "\t", na = ".")
```

## `0_Rscripts/5_var_density.r`

Raw source: lines 320–601 of the compilation.

```r
# ============================================================
# Sliding window AF(R) - AF(S) + densité de variants (plotly)
# - Fenêtre = 500000 bp
# - Pas = 100000 bp
# - Un panel par chromosome (colonne 1 = chromosome)
# - Aire sous la courbe colorée par densité (viridis)
# - Hover sur les points de la courbe : position + valeurs
# - Remplissage en trapèzes non chevauchants sous la courbe
# ============================================================

# Packages
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(plotly)
  library(htmlwidgets)
})

# ---- Paramètres ----
tsv_file <- "./4_join/var.tsv"   # <- adapte le chemin
window_size <- 500000
step_size   <- 10000

# ---- Lecture ----
dt <- fread(tsv_file, sep = "\t", header = TRUE, na.strings = c("", "NA"))

# Vérifs minimales
required_cols <- c("chromosome", "variant_pos")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Colonnes manquantes : ", paste(missing_cols, collapse = ", "))
}

# Colonnes AF
s_cols <- grep("^S_", names(dt), value = TRUE)
r_cols <- grep("^R_", names(dt), value = TRUE)

if (length(s_cols) == 0 || length(r_cols) == 0) {
  stop("Impossible de trouver les colonnes AF S_... et/ou R_...")
}

# ---- Nettoyage / typage ----
dt[, variant_pos := suppressWarnings(as.numeric(variant_pos))]

# Convertit les AF en numérique
for (cc in c(s_cols, r_cols)) {
  dt[, (cc) := suppressWarnings(as.numeric(get(cc)))]
}

# Garde les lignes avec position valide
dt <- dt[!is.na(variant_pos)]

# ---- Moyenne AF par variant ----
# rowMeans(na.rm=TRUE) renvoie NaN si toute la ligne est NA -> on remet NA
dt[, mean_AF_S := rowMeans(.SD, na.rm = TRUE), .SDcols = s_cols]
dt[, mean_AF_R := rowMeans(.SD, na.rm = TRUE), .SDcols = r_cols]

dt[is.nan(mean_AF_S), mean_AF_S := NA_real_]
dt[is.nan(mean_AF_R), mean_AF_R := NA_real_]

# Différence demandée : moyenne(R) - moyenne(S)
dt[, delta_AF := mean_AF_R - mean_AF_S]

# On garde les variants informatifs
dtv <- dt[!is.na(delta_AF), .(chromosome, variant_pos, delta_AF)]

if (nrow(dtv) == 0) {
  stop("Aucun variant avec delta_AF calculable.")
}

# ---- Fonction sliding window par chromosome ----
make_windows_for_chr <- function(dchr, chr_name, window_size = 5e5, step_size = 1e5) {
  min_pos <- floor(min(dchr$variant_pos, na.rm = TRUE))
  max_pos <- ceiling(max(dchr$variant_pos, na.rm = TRUE))
  
  # On démarre à 1 (ou au min) ; ici on prend 1 pour une échelle génomique plus lisible
  starts <- seq(1, max_pos, by = step_size)
  
  wins <- data.table(
    chromosome = chr_name,
    win_start = starts,
    win_end   = starts + window_size - 1L
  )
  wins[, win_id := .I]
  
  # Prépare les intervals pour foverlaps
  vars_int <- copy(dchr)
  vars_int[, `:=`(start = variant_pos, end = variant_pos)]
  setkey(vars_int, start, end)
  
  wins_int <- copy(wins)
  setkey(wins_int, win_start, win_end)
  
  # Overlap variants x windows
  ov <- foverlaps(
    x = vars_int,
    y = wins_int,
    by.x = c("start", "end"),
    by.y = c("win_start", "win_end"),
    type = "within",
    nomatch = 0L
  )
  
  # Agrégation par fenêtre
  agg <- ov[, .(
    n_variants = .N,
    mean_delta = mean(delta_AF, na.rm = TRUE)
  ), by = .(chromosome, win_id, win_start, win_end)]
  
  # Ajoute les fenêtres vides
  out <- merge(
    wins[, .(chromosome, win_id, win_start, win_end)],
    agg,
    by = c("chromosome", "win_id", "win_start", "win_end"),
    all.x = TRUE
  )
  
  out[is.na(n_variants), n_variants := 0L]
  out[, win_mid := (win_start + win_end) / 2]
  
  out[]
}

# ---- Calcul sur tous les chromosomes ----
windows_list <- lapply(split(dtv, by = "chromosome", keep.by = TRUE), function(x) {
  make_windows_for_chr(x, chr_name = x$chromosome[1], window_size = window_size, step_size = step_size)
})

win_dt <- rbindlist(windows_list, use.names = TRUE, fill = TRUE)

# ---- Préparation de l'aire colorée sous la courbe (baseline par chromosome) ----
chr_ranges <- win_dt[, .(
  y_min = suppressWarnings(min(mean_delta, na.rm = TRUE)),
  y_max = suppressWarnings(max(mean_delta, na.rm = TRUE))
), by = chromosome]

# Cas pathologiques (tout NA)
chr_ranges[!is.finite(y_min) | !is.finite(y_max), `:=`(y_min = -0.1, y_max = 0.1)]
chr_ranges[, y_span := y_max - y_min]
chr_ranges[y_span == 0 | !is.finite(y_span), y_span := 0.2]

# Ligne de base pour l'aire (sous la courbe)
chr_ranges[, baseline := y_min - 0.15 * y_span]

win_dt <- merge(
  win_dt,
  chr_ranges[, .(chromosome, baseline)],
  by = "chromosome",
  all.x = TRUE
)

# Texte de hover (sur les points de la courbe)
win_dt[, hover_txt := paste0(
  "Chromosome: ", chromosome,
  "<br>Window start: ", format(win_start, big.mark = " "),
  "<br>Window end: ", format(win_end, big.mark = " "),
  "<br>Window mid: ", format(round(win_mid), big.mark = " "),
  "<br>Mean AF(R)-AF(S): ", ifelse(is.na(mean_delta), "NA", sprintf("%.4f", mean_delta)),
  "<br>Variant density: ", n_variants
)]

# ---- Construction de trapèzes non chevauchants sous la courbe ----
# Chaque trapèze relie 2 points consécutifs de la courbe :
# (x_i, baseline) -> (x_i, y_i) -> (x_{i+1}, y_{i+1}) -> (x_{i+1}, baseline)
build_trapezoids <- function(dchr) {
  d <- copy(dchr[!is.na(mean_delta)])
  setorder(d, win_mid)
  
  if (nrow(d) < 2) return(data.table())
  
  # Segment i -> i+1
  seg <- d[, .(
    chromosome,
    x_left   = win_mid,
    y_left   = mean_delta,
    dens_left = n_variants,
    baseline
  )]
  seg[, `:=`(
    x_right   = shift(x_left, type = "lead"),
    y_right   = shift(y_left, type = "lead"),
    dens_right = shift(dens_left, type = "lead")
  )]
  seg <- seg[!is.na(x_right)]
  
  if (nrow(seg) == 0) return(data.table())
  
  # Couleur du trapèze = moyenne des densités des 2 points
  seg[, fill_density := (dens_left + dens_right) / 2]
  seg[, trap_id := paste0(chromosome, "__", .I)]
  
  # 4 sommets par trapèze
  poly <- rbindlist(list(
    seg[, .(chromosome, trap_id, ord = 1L, x = x_left,  y = baseline,     fill_density)],
    seg[, .(chromosome, trap_id, ord = 2L, x = x_left,  y = y_left,       fill_density)],
    seg[, .(chromosome, trap_id, ord = 3L, x = x_right, y = y_right,      fill_density)],
    seg[, .(chromosome, trap_id, ord = 4L, x = x_right, y = baseline,     fill_density)]
  ), use.names = TRUE)
  
  setorder(poly, chromosome, trap_id, ord)
  poly[]
}

poly_list <- lapply(split(win_dt, by = "chromosome", keep.by = TRUE), build_trapezoids)
poly_dt <- rbindlist(poly_list, use.names = TRUE, fill = TRUE)

# ---- Plot ggplot (base) ----
p <- ggplot() +
  # Trapèzes non chevauchants sous la courbe colorés par densité
  geom_polygon(
    data = poly_dt,
    aes(
      x = x,
      y = y,
      group = trap_id,
      fill = fill_density
    ),
    alpha = 0.9,
    colour = NA
  ) +
  # Courbe delta AF
  geom_line(
    data = win_dt[!is.na(mean_delta)],
    aes(x = win_mid, y = mean_delta, group = chromosome),
    linewidth = 0.7,
    color = "black"
  ) +
  # Points pour hover plotly
  geom_point(
    data = win_dt[!is.na(mean_delta)],
    aes(x = win_mid, y = mean_delta, text = hover_txt),
    size = 0.8,
    alpha = 0.01   # quasi invisible, mais survolable
  ) +
  # Ligne y = 0
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey40"
  ) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 1) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Variant density\n(count/window)"
  ) +
  labs(
    title = "Sliding window: mean AF(R) - mean AF(S)",
    subtitle = paste0("Window = ", format(window_size, big.mark = " "), " bp | Step = ", format(step_size, big.mark = " "), " bp"),
    x = "Genomic position (bp)",
    y = "Mean AF(R) - Mean AF(S)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# ---- Conversion en plotly (plus haut, moins large) ----
n_chr <- length(unique(win_dt$chromosome))

pp <- ggplotly(
  p,
  tooltip = "text",
  width = 900,                    # moins large
  height = max(700, 280 * n_chr)  # plus haut
) %>%
  layout(
    hovermode = "closest"
  )

# Affichage interactif (utile en session R/RStudio)
print(pp)

# ---- Sauvegarde HTML interactive (IMPORTANT avec Rscript) ----
htmlwidgets::saveWidget(
  widget = pp,
  file = "./5_var_density/sliding_AF_R_minus_S_by_chr_plotly.html",
  selfcontained = TRUE
)
```

## `0_Rscripts/6_plotly_var.r`

Raw source: lines 607–991 of the compilation.

```r
#!/usr/bin/env Rscript

## ==== Packages ==== ##
library(readr)
library(dplyr)
library(stringr)
library(plotly)
library(htmlwidgets)
library(rmarkdown)

## ==== 1. Lecture du fichier ==== ##
df <- read_tsv("./4_join/var.tsv", show_col_types = FALSE)

## ==== 1bis. Harmonisation minimale de noms de colonnes (si besoin) ==== ##
if (!"variant_ref" %in% names(df) && "ref" %in% names(df)) {
  df <- df %>% mutate(variant_ref = ref)
}
if (!"variant_alt" %in% names(df) && "alt" %in% names(df)) {
  df <- df %>% mutate(variant_alt = alt)
}
if (!"variant_pos" %in% names(df) && "pos" %in% names(df)) {
  df <- df %>% mutate(variant_pos = pos)
}

## ==== 1ter. Colonnes positions numériques (si présentes) ==== ##
for (cn in intersect(c("variant_pos", "gene_start", "gene_end"), names(df))) {
  df[[cn]] <- suppressWarnings(as.numeric(df[[cn]]))
}

## ==== 2. Détection des colonnes d'AF ==== ##
af_cols   <- grep("^[RS]", names(df), value = TRUE)
r_samples <- grep("^R", names(df), value = TRUE)
s_samples <- grep("^S", names(df), value = TRUE)

## ==== 2bis. Colonnes d'impact (compteurs) ==== ##
impact_count_cols <- intersect(c("HIGH", "MODERATE", "MODIFIER", "LOW"), names(df))

## ==== 2ter. Colonnes potentielles d'effets prédits (pour hover) ==== ##
effect_candidate_cols <- intersect(
  c(
    "predicted_effect", "predicted_effects", "effect", "effects",
    "annotation", "annotations", "consequence", "consequences",
    "impact", "impacts", "ANN", "CSQ", "snpEff", "snpeff", "snpeff_effect"
  ),
  names(df)
)

## ==== 3. Conversion des colonnes numériques (AF + impacts) ==== ##
df <- df %>%
  mutate(
    across(all_of(af_cols), ~ {
      x <- as.character(.x)
      x[x == "."] <- NA_character_
      x <- str_replace(x, ",", ".")
      suppressWarnings(as.numeric(x))
    }),
    across(all_of(impact_count_cols), ~ suppressWarnings(as.numeric(.x)))
  )

## ==== 3bis. Helpers texte (effets prédits / AF par sample) ==== ##
# Effets prédits concaténés (si colonnes disponibles)
if (length(effect_candidate_cols) > 0) {
  eff_df <- df %>% select(all_of(effect_candidate_cols))
  df$predicted_effects_info <- apply(eff_df, 1, function(v) {
    v <- as.character(v)
    v <- v[!is.na(v) & v != "" & v != "."]
    v <- unique(v)
    if (length(v) == 0) "NA" else paste(v, collapse = " | ")
  })
} else {
  df$predicted_effects_info <- "NA"
}

# AF par sample pour hover
if (length(af_cols) > 0) {
  af_df <- df %>% select(all_of(af_cols))
  df$af_info <- apply(af_df, 1, function(v) {
    vals <- ifelse(is.na(v), "NA", sprintf("%.3f", as.numeric(v)))
    paste0(paste0(af_cols, " = ", vals), collapse = "<br>")
  })
} else {
  df$af_info <- "No AF columns detected"
}

## ==== 4. Préparation au niveau variant (une bulle = un variant) ==== ##
# Identifiant variant robuste
df <- df %>%
  mutate(
    variant_id = paste(
      coalesce(as.character(chromosome), "NA"),
      coalesce(as.character(variant_pos), "NA"),
      coalesce(as.character(variant_ref), "NA"),
      coalesce(as.character(variant_alt), "NA"),
      sep = ":"
    )
  )

# Impact "meilleur" par variant (même logique qu'avant, mais ligne par ligne)
if (all(c("HIGH", "MODERATE", "MODIFIER", "LOW") %in% names(df))) {
  df <- df %>%
    mutate(
      best_impact = case_when(
        !is.na(HIGH)     & HIGH > 0     ~ "HIGH",
        !is.na(MODERATE) & MODERATE > 0 ~ "MODERATE",
        !is.na(MODIFIER) & MODIFIER > 0 ~ "MODIFIER",
        !is.na(LOW)      & LOW > 0      ~ "LOW",
        TRUE                          ~ "NONE"
      )
    )
} else {
  # fallback: si pas de compteurs, on reste sur NONE
  df <- df %>% mutate(best_impact = "NONE")
}

# Moyennes AF R / S et diff
if (length(r_samples) > 0) {
  df$mean_R <- rowMeans(df[, r_samples, drop = FALSE], na.rm = TRUE)
  df$mean_R[!is.finite(df$mean_R)] <- NA_real_
} else {
  df$mean_R <- NA_real_
}

if (length(s_samples) > 0) {
  df$mean_S <- rowMeans(df[, s_samples, drop = FALSE], na.rm = TRUE)
  df$mean_S[!is.finite(df$mean_S)] <- NA_real_
} else {
  df$mean_S <- NA_real_
}

df <- df %>%
  mutate(
    diff_AF = mean_R - mean_S
  )

# Infos impacts (compteurs) pour hover
for (cn in c("HIGH", "MODERATE", "MODIFIER", "LOW")) {
  if (!cn %in% names(df)) df[[cn]] <- NA_real_
}

# Colonnes exoniques potentielles (si présentes)
exon_start_col <- names(df)[tolower(names(df)) %in% c("exon_start", "exonstart")]
exon_end_col   <- names(df)[tolower(names(df)) %in% c("exon_end", "exonend")]

# Hover texte complet par variant
df <- df %>%
  mutate(
    gene_label = ifelse(!is.na(gene) & gene != "", as.character(gene), "Intergenic / NA"),
    region_label = if ("region_code" %in% names(df)) as.character(region_code) else "NA",
    gene_start_label = if ("gene_start" %in% names(df)) ifelse(is.na(gene_start), "NA", as.character(gene_start)) else "NA",
    gene_end_label   = if ("gene_end" %in% names(df)) ifelse(is.na(gene_end), "NA", as.character(gene_end)) else "NA",
    exon_start_label = if (length(exon_start_col) == 1) ifelse(is.na(.data[[exon_start_col]]), "NA", as.character(.data[[exon_start_col]])) else "NA",
    exon_end_label   = if (length(exon_end_col) == 1) ifelse(is.na(.data[[exon_end_col]]), "NA", as.character(.data[[exon_end_col]])) else "NA",
    hover_text = paste0(
      "Chromosome: ", coalesce(as.character(chromosome), "NA"),
      "<br>Position: ", coalesce(as.character(variant_pos), "NA"),
      "<br>Variant: ", coalesce(as.character(variant_ref), "NA"), ">", coalesce(as.character(variant_alt), "NA"),
      "<br>Region code: ", coalesce(region_label, "NA"),
      "<br>Gene: ", gene_label,
      "<br>Gene start-end: ", gene_start_label, " - ", gene_end_label,
      "<br>Exon start-end: ", exon_start_label, " - ", exon_end_label,
      "<br>Predicted effects: ", predicted_effects_info,
      "<br>Impact counts: H=", coalesce(as.character(HIGH), "NA"),
      ", M=", coalesce(as.character(MODERATE), "NA"),
      ", Md=", coalesce(as.character(MODIFIER), "NA"),
      ", L=", coalesce(as.character(LOW), "NA"),
      "<br>Best impact: ", best_impact,
      "<br>mean_R: ", ifelse(is.na(mean_R), "NA", sprintf("%.3f", mean_R)),
      "<br>mean_S: ", ifelse(is.na(mean_S), "NA", sprintf("%.3f", mean_S)),
      "<br>diff_AF (R - S): ", ifelse(is.na(diff_AF), "NA", sprintf("%.3f", diff_AF)),
      "<br><br>AF by sample:<br>", af_info
    )
  ) %>%
  mutate(
    best_impact = factor(
      best_impact,
      levels = c("HIGH", "MODERATE", "MODIFIER", "LOW", "NONE")
    )
  )

## ==== 4bis. Export table : top 10% des variants (|diff_AF|) ==== ##
top10_threshold <- quantile(abs(df$diff_AF), probs = 0.90, na.rm = TRUE)

top10_variants <- df %>%
  filter(!is.na(diff_AF)) %>%
  mutate(abs_diff_AF = abs(diff_AF)) %>%
  filter(abs_diff_AF >= top10_threshold) %>%
  arrange(desc(abs_diff_AF))

write_tsv(top10_variants, "./6_plotly_var/top10pct_variants_by_absDiffAF.tsv")

## ==== 4ter. Préparation des plages gènes / exons pour arrière-plan ==== ##
# Gènes (rectangles verticaux légers)
gene_ranges <- tibble()
if (all(c("chromosome", "gene", "gene_start", "gene_end") %in% names(df))) {
  gene_ranges <- df %>%
    filter(
      !is.na(chromosome),
      !is.na(gene),
      gene != "",
      !is.na(gene_start),
      !is.na(gene_end)
    ) %>%
    distinct(chromosome, gene, gene_start, gene_end) %>%
    filter(gene_end >= gene_start)
}

# Exons (si colonnes exon_start/exon_end existent)
exon_ranges <- tibble()
if (length(exon_start_col) == 1 && length(exon_end_col) == 1) {
  exon_ranges <- df %>%
    filter(
      !is.na(chromosome),
      !is.na(.data[[exon_start_col]]),
      !is.na(.data[[exon_end_col]])
    ) %>%
    transmute(
      chromosome = chromosome,
      exon_start = suppressWarnings(as.numeric(.data[[exon_start_col]])),
      exon_end   = suppressWarnings(as.numeric(.data[[exon_end_col]]))
    ) %>%
    filter(!is.na(exon_start), !is.na(exon_end), exon_end >= exon_start) %>%
    distinct()
} else if ("region_code" %in% names(df) && "variant_pos" %in% names(df)) {
  # Fallback : si pas de bornes exoniques, utiliser les variants annotés exon ("e") comme traits fins
  exon_ranges <- df %>%
    filter(!is.na(chromosome), region_code == "e", !is.na(variant_pos)) %>%
    transmute(
      chromosome = chromosome,
      exon_start = variant_pos,
      exon_end   = variant_pos
    ) %>%
    distinct()
}

## ==== 5. Génération d'un plotly par chromosome ==== ##
plots_by_chr <- lapply(
  split(df, df$chromosome),
  function(dat_chr) {

    layers_order <- c("NONE", "LOW", "MODIFIER", "MODERATE", "HIGH")
    chr_name <- unique(dat_chr$chromosome)[1]

    # bornes x (variants présents)
    x_min <- suppressWarnings(min(dat_chr$variant_pos, na.rm = TRUE))
    x_max <- suppressWarnings(max(dat_chr$variant_pos, na.rm = TRUE))
    if (!is.finite(x_min) || !is.finite(x_max)) {
      x_min <- 0
      x_max <- 1
    }

    # shapes de fond : gènes puis exons
    shapes_bg <- list()

    # Gènes (fond très léger)
    if (nrow(gene_ranges) > 0) {
      g_chr <- gene_ranges %>% filter(chromosome == chr_name)
      if (nrow(g_chr) > 0) {
        for (i in seq_len(nrow(g_chr))) {
          shapes_bg[[length(shapes_bg) + 1]] <- list(
            type = "rect",
            xref = "x", yref = "paper",
            x0 = g_chr$gene_start[i], x1 = g_chr$gene_end[i],
            y0 = 0, y1 = 1,
            fillcolor = "rgba(70,130,180,0.08)",  # bleu léger
            line = list(width = 0),
            layer = "below"
          )
        }
      }
    }

    # Exons (fond encore visible mais léger)
    if (nrow(exon_ranges) > 0) {
      e_chr <- exon_ranges %>% filter(chromosome == chr_name)
      if (nrow(e_chr) > 0) {
        for (i in seq_len(nrow(e_chr))) {
          shapes_bg[[length(shapes_bg) + 1]] <- list(
            type = "rect",
            xref = "x", yref = "paper",
            x0 = e_chr$exon_start[i], x1 = e_chr$exon_end[i],
            y0 = 0, y1 = 1,
            fillcolor = "rgba(34,139,34,0.12)",   # vert léger
            line = list(width = 0),
            layer = "below"
          )
        }
      }
    }

    p <- plot_ly()

    # Traces par impact (une bulle = un variant, taille fixe)
    for (imp in layers_order) {
      dat_layer <- dat_chr %>% filter(best_impact == imp)
      if (nrow(dat_layer) == 0) next

      col_rgba <- switch(
        imp,
        "HIGH"     = "rgba(187,11,11,0.95)",
        "MODERATE" = "rgba(255,127,0,0.95)",
        "MODIFIER" = "rgba(199,207,0,0.95)",
        "LOW"      = "rgba(20,148,20,0.95)",
        "NONE"     = "rgba(128,128,128,0.90)",
        "rgba(128,128,128,0.90)"
      )

      p <- p %>%
        add_trace(
          data = dat_layer,
          x = ~variant_pos,
          y = ~diff_AF,
          type = "scatter",
          mode = "markers",
          text = ~hover_text,
          hoverinfo = "text",
          name = imp,
          legendgroup = imp,
          marker = list(
            size = 7,
            color = col_rgba,
            line = list(width = 0)
          ),
          showlegend = TRUE
        )
    }

    # Labels de gènes (optionnels, masqués par défaut), positionnés au milieu des gènes
    if (nrow(gene_ranges) > 0) {
      gene_lab_chr <- gene_ranges %>%
        filter(chromosome == chr_name) %>%
        mutate(gene_mid = (gene_start + gene_end) / 2)

      if (nrow(gene_lab_chr) > 0) {
        p <- p %>%
          add_trace(
            data = gene_lab_chr,
            x = ~gene_mid,
            y = 0,
            type = "scatter",
            mode = "text",
            text = ~gene,
            textposition = "middle center",
            textfont = list(color = "rgba(0,0,0,0.25)", size = 9),
            hoverinfo = "skip",
            name = "Gene labels",
            showlegend = TRUE,
            visible = "legendonly"
          )
      }
    }

    p %>%
      layout(
        title = chr_name,
        xaxis = list(title = "Variant position (bp)", range = c(x_min, x_max)),
        yaxis = list(title = "mean(AF_R) - mean(AF_S)"),
        legend = list(title = list(text = "Impact / Labels")),
        shapes = shapes_bg
      )
  }
)

## ==== 6. Sauvegarde des plots dans des fichiers HTML ==== ##
if (!rmarkdown::pandoc_available()) {
  stop(
    "Pandoc introuvable : impossible de produire des HTML self-contained.\n",
    "Installe Pandoc (ex: 'sudo apt-get install pandoc' sur Debian/Ubuntu,\n",
    "ou via conda: 'conda install -c conda-forge pandoc'), puis relance."
  )
}

file_names <- character(0)

for (chr in names(plots_by_chr)) {
  p <- plots_by_chr[[chr]]
  safe_chr <- gsub("[^A-Za-z0-9_-]", "_", chr)
  file_name <- paste0("./6_plotly_var/", safe_chr, ".html")
  saveWidget(p, file = file_name, selfcontained = TRUE)
  file_names <- c(file_names, file_name)
}

cat("Fichiers générés :\n")
cat(paste0(" - ", file_names), sep = "\n")
cat("\nTable exportée :\n")
cat(" - top10pct_variants_by_absDiffAF.tsv\n")
```

## `0_Rscripts/7_ggplot_gene.r`

Raw source: lines 997–1313 of the compilation.

```r
#!/usr/bin/env Rscript

## ==== Packages ====
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)

## ==== 0. Récupération des bornes optionnelles (-l / -r / -b / -n / -f) ====
args <- commandArgs(trailingOnly = TRUE)

get_opt_num <- function(flag) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) {
    suppressWarnings(as.numeric(args[i + 1]))
  } else {
    NA_real_
  }
}

has_flag <- function(flag) {
  any(args == flag)
}

l <- get_opt_num("-l")   # borne gauche (bp)
r <- get_opt_num("-r")   # borne droite (bp)
b <- get_opt_num("-b")   # borne basse (diff_AF)
n <- has_flag("-n")      # ajoute labels gènes
f <- has_flag("-f")      # force y-axis de 0 à 1

## ==== 1. Lecture du fichier (comme 7_plotly_gene.r) ====
df <- read_tsv("./4_join/var.tsv", show_col_types = FALSE)

## ==== 1bis. Garder uniquement les variants dans des gènes (g ou e) ====
if ("region_code" %in% names(df)) {
  df <- df %>% filter(region_code %in% c("g", "e"))
}

## ==== 1ter. Harmonisation minimale de noms de colonnes (si besoin) ====
if (!"variant_ref" %in% names(df) && "ref" %in% names(df)) {
  df <- df %>% mutate(variant_ref = ref)
}
if (!"variant_alt" %in% names(df) && "alt" %in% names(df)) {
  df <- df %>% mutate(variant_alt = alt)
}

## ==== 2. Détection des colonnes d'AF ====
# colonnes AF = commencent par R ou S (RARA, RBET, SCHI, SLUC, ...)
af_cols <- grep("^[RS]", names(df), value = TRUE)
r_samples <- grep("^R", names(df), value = TRUE)
s_samples <- grep("^S", names(df), value = TRUE)

## ==== 2bis. Colonnes d'impact (compteurs) ====
impact_count_cols <- intersect(c("HIGH", "MODERATE", "MODIFIER", "LOW"), names(df))
has_impact_counts <- length(impact_count_cols) > 0
has_snpeff_impact <- "snpeff_impact" %in% names(df)
has_snpeff_effect <- "snpeff_effect" %in% names(df)

## ==== 3. Conversion des colonnes numériques (AF + impacts) ====
df <- df %>%
  mutate(
    across(all_of(af_cols), ~ {
      x <- as.character(.x)   # <- clé: uniformiser le type
      x[x == "."] <- NA_character_
      x <- str_replace(x, ",", ".")
      as.numeric(x)
    }),
    across(all_of(impact_count_cols), ~ suppressWarnings(as.numeric(.x)))
  )

## ==== 4. Agrégation au niveau gène ====
gene_summary <- df %>%
  group_by(chromosome, gene, gene_start, gene_end) %>%
  summarise(
    ## position milieu du gène
    gene_mid = (first(gene_start) + first(gene_end)) / 2,

    ## longueur du gène (bp) - robuste même si start/end inversés
    gene_len_bp = abs(first(gene_end) - first(gene_start)) + 1,

    ## nombre de variants uniques pos/ref/alt
    n_variants = n_distinct(paste(variant_pos, variant_ref, variant_alt, sep = ":")),

    ## AF moyenne par sample sur tous les variants du gène
    across(all_of(af_cols), ~ mean(.x, na.rm = TRUE)),

    ## meilleur impact : comme 7_plotly_gene.r si compteurs présents, sinon fallback SnpEff si dispo
    best_impact = dplyr::case_when(
      has_impact_counts && ("HIGH" %in% impact_count_cols) && any(HIGH > 0, na.rm = TRUE) ~ "HIGH",
      has_impact_counts && ("MODERATE" %in% impact_count_cols) && any(MODERATE > 0, na.rm = TRUE) ~ "MODERATE",
      has_impact_counts && ("MODIFIER" %in% impact_count_cols) && any(MODIFIER > 0, na.rm = TRUE) ~ "MODIFIER",
      has_impact_counts && ("LOW" %in% impact_count_cols) && any(LOW > 0, na.rm = TRUE) ~ "LOW",
      (!has_impact_counts) && has_snpeff_impact && any(snpeff_impact == "HIGH", na.rm = TRUE) ~ "HIGH",
      (!has_impact_counts) && has_snpeff_impact && any(snpeff_impact == "MODERATE", na.rm = TRUE) ~ "MODERATE",
      (!has_impact_counts) && has_snpeff_impact && any(snpeff_impact == "MODIFIER", na.rm = TRUE) ~ "MODIFIER",
      (!has_impact_counts) && has_snpeff_impact && any(snpeff_impact == "LOW", na.rm = TRUE) ~ "LOW",
      TRUE ~ "NONE"
    ),

    ## concat de toutes les infos variants pos/ref/alt + impacts (comme plotly si compteurs présents)
    variant_info = if (has_impact_counts) {
      paste(
        unique(paste0(
          variant_pos, " ", variant_ref, ">", variant_alt,
          " [H=", coalesce(as.character(if ("HIGH" %in% impact_count_cols) HIGH else NA), "NA"),
          ", M=", coalesce(as.character(if ("MODERATE" %in% impact_count_cols) MODERATE else NA), "NA"),
          ", Md=", coalesce(as.character(if ("MODIFIER" %in% impact_count_cols) MODIFIER else NA), "NA"),
          ", L=", coalesce(as.character(if ("LOW" %in% impact_count_cols) LOW else NA), "NA"),
          "]"
        )),
        collapse = "; "
      )
    } else if (has_snpeff_effect && has_snpeff_impact) {
      paste(
        unique(paste0(
          variant_pos, " ", variant_ref, ">", variant_alt,
          " [", snpeff_effect, ":", snpeff_impact, "]"
        )),
        collapse = "; "
      )
    } else {
      paste(
        unique(paste0(
          variant_pos, " ", variant_ref, ">", variant_alt
        )),
        collapse = "; "
      )
    },
    .groups = "drop"
  ) %>%
  ## calculs par ligne (par gène)
  rowwise() %>%
  mutate(
    ## moyennes AF R vs S
    mean_R = mean(c_across(all_of(r_samples)), na.rm = TRUE),
    mean_S = mean(c_across(all_of(s_samples)), na.rm = TRUE),
    diff_AF = mean_R - mean_S,

    ## densité de mutations normalisée par taille de gène
    ## (variants par kb = n_variants / (gene_len_bp/1000))
    variants_per_kb = n_variants * 1000 / gene_len_bp,

    ## texte AF par sample (conservé, même si ggplot n'a pas de hover natif)
    af_info = paste(
      paste0(
        af_cols, " = ", sprintf("%.3f", c_across(all_of(af_cols)))
      ),
      collapse = "\n"
    ),

    ## texte complet (conservé pour exports/inspection)
    hover_text = paste0(
      "Chromosome: ", chromosome,
      "\nGene: ", gene,
      "\nStart: ", gene_start,
      "\nEnd: ", gene_end,
      "\nMid pos: ", round(gene_mid),
      "\nGene length (bp): ", gene_len_bp,
      "\n#variants: ", n_variants,
      "\nVariants per kb: ", sprintf("%.3f", variants_per_kb),
      "\nBest impact: ", best_impact,
      "\nmean_R: ", sprintf("%.3f", mean_R),
      "\nmean_S: ", sprintf("%.3f", mean_S),
      "\ndiff_AF (R - S): ", sprintf("%.3f", diff_AF),
      "\n\nAF by sample:\n", af_info,
      "\n\nVariants:\n", variant_info
    )
  ) %>%
  ungroup() %>%
  mutate(
    best_impact = factor(
      best_impact,
      levels = c("HIGH", "MODERATE", "LOW", "MODIFIER", "NONE")
    )
  )

## ==== 4bis. Export table : top 10% des gènes (|diff_AF|) ====
# Définition "top 10%" = les gènes dont |diff_AF| est dans le décile supérieur
# (si tu veux plutôt "diff_AF positif seulement", dis-le et je te le ajuste)
top10_threshold <- quantile(abs(gene_summary$diff_AF), probs = 0.90, na.rm = TRUE)

top10_genes <- gene_summary %>%
  filter(!is.na(diff_AF)) %>%
  mutate(abs_diff_AF = abs(diff_AF)) %>%
  filter(abs_diff_AF >= top10_threshold) %>%
  arrange(desc(abs_diff_AF))

# Export TSV (inclut maintenant gene_len_bp, variants_per_kb)
write_tsv(top10_genes, "./7_ggplot_gene/top10pct_genes_by_absDiffAF.tsv")

## ==== 5. Génération d'un ggplot par chromosome ====
plots_by_chr <- lapply(
  split(gene_summary, gene_summary$chromosome),
  function(dat_chr) {

    ## bornes par chromosome (bp)
    chr_max <- suppressWarnings(max(dat_chr$gene_mid, na.rm = TRUE))
    left_bp <- if (!is.na(l)) max(0, l) else 0
    right_bp <- if (!is.na(r)) min(r, chr_max) else chr_max

    ## borne basse (y) si demandée
    bottom_y <- if (!is.na(b)) b else NA_real_

    ## on filtre à l'intervalle si demandé (X)
    dat_chr <- dat_chr %>%
      filter(gene_mid >= left_bp, gene_mid <= right_bp)

    ## sous-dataframe strictement "dans la fenêtre" pour les étiquettes
    ## (évite d'étiqueter des points hors zone, surtout quand ylim est rogné)
    dat_lbl <- dat_chr
    if (!is.na(bottom_y)) {
      dat_lbl <- dat_lbl %>% filter(!is.na(diff_AF), diff_AF >= bottom_y)
    }
    if (f) {
      dat_lbl <- dat_lbl %>% filter(!is.na(diff_AF), diff_AF >= 0, diff_AF <= 1)
    }

    ## ordre d’affichage: ce qui est dessiné en dernier est au-dessus
    ## (LOW tout en dessous, puis MODIFIER, puis MODERATE, puis HIGH au-dessus)
    dat_chr <- dat_chr %>%
      mutate(
        .draw_order = dplyr::case_when(
          best_impact == "NONE" ~ 1,
          best_impact == "LOW" ~ 2,
          best_impact == "MODIFIER" ~ 3,
          best_impact == "MODERATE" ~ 4,
          best_impact == "HIGH" ~ 5,
          TRUE ~ 0
        )
      ) %>%
      arrange(.draw_order)

    ggplot(
      dat_chr,
      aes(
        x = gene_mid / 1e6,
        y = diff_AF,
        size = variants_per_kb,   # <-- taille = variants/kb (sans log)
        color = best_impact
      )
    ) +
      geom_point(shape = 16, alpha = 0.45) +
      coord_cartesian(
        xlim = c(left_bp, right_bp) / 1e6,
        ylim = if (f) c(0, 1) else if (!is.na(bottom_y)) c(bottom_y, NA) else NULL
      ) +
      scale_color_manual(
        values = c(
          HIGH = "#440154",       # CC0000
          MODERATE = "#440154",   # FF6600
          LOW = "#440154",        # 339966
          MODIFIER = "#440154",   # CCCC00
          NONE = "#440154"        # grey40
        ),
        drop = FALSE
      ) +
      scale_size_continuous(range = c(0.1, 4.5), breaks = waiver()) +
      labs(
        title = paste0("Chromosome ", unique(dat_chr$chromosome)),
        x = "Genomic position (Mb)",
        y = "µAF_R-µAF_S",
        color = "Best impact",
        size = "Variants / kb"
      ) +
      theme_bw(base_size = 24) +
      theme(
        plot.title = element_text(face = "bold", size = 24, color = "black", hjust = 0.5),
        axis.title.x = element_text(size = 24, color = "black"),
        axis.title.y = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 24, color = "black"),
        axis.text.y = element_text(size = 24, color = "black"),
        legend.title = element_text(face = "bold", size = 24, color = "black"),
        legend.text = element_text(size = 24, color = "black"),
        panel.grid.minor = element_blank()
      ) +
      {
        if (n) geom_text_repel(
          data = dat_lbl,
          aes(label = gene),
          size = 3,                # petit
          min.segment.length = 0, # segment toujours dessiné si possible
          segment.size = 0.3,
          box.padding = 0.35,
          point.padding = 0.25,
          max.overlaps = Inf
        )
      }
  }
)

## ==== 6. Sauvegarde des plots dans des fichiers PDF ====
file_names <- character(0)

for (chr in names(plots_by_chr)) {
  p <- plots_by_chr[[chr]]

  # Nom de fichier safe (remplacer caractères bizarres par "_")
  safe_chr <- gsub("[^A-Za-z0-9_-]", "_", chr)
  file_name <- paste0("./7_ggplot_gene/ggplot_", safe_chr, ".pdf")

  ggsave(
    filename = file_name,
    plot = p,
    width = 20,
    height = 5,
    units = "in",
    bg = "white"
  )

  file_names <- c(file_names, file_name)
}

cat("Fichiers générés :\n")
cat(paste0(" - ", file_names), sep = "\n")
cat("\nTable exportée :\n")
cat(" - top10pct_genes_by_absDiffAF.tsv\n")
```

## `0_Rscripts/7_plotly_gene.r`

Raw source: lines 1319–1742 of the compilation.

```r
#!/usr/bin/env Rscript

## ==== Packages ==== ##
library(readr)
library(dplyr)
library(stringr)
library(plotly)
library(htmlwidgets)
library(rmarkdown)

## ==== 1. Lecture du fichier ==== ##
df <- read_tsv("./4_join/var.tsv", show_col_types = FALSE)

## ==== 1bis. Garder uniquement les variants dans des gènes (g ou e) ==== ##
df <- df %>%
  filter(region_code %in% c("g", "e"))

## ==== 1ter. Harmonisation minimale de noms de colonnes (si besoin) ==== ##
if (!"variant_ref" %in% names(df) && "ref" %in% names(df)) {
  df <- df %>% mutate(variant_ref = ref)
}
if (!"variant_alt" %in% names(df) && "alt" %in% names(df)) {
  df <- df %>% mutate(variant_alt = alt)
}

## ==== 1quater. Fonctions rapides pour INFO/ANN ==== ##

clean_str <- function(x) {
  if (length(x) == 0 || is.na(x)) return(NA_character_)
  x <- trimws(as.character(x))
  if (x %in% c("", ".", "NA", "na")) return(NA_character_)
  x
}

clean_hgvsp <- function(x) {
  x <- clean_str(x)
  if (is.na(x)) return(NA_character_)
  x <- sub("^p\\.", "", x)
  x <- clean_str(x)
  x
}

effect_label_fast <- function(annotation, hgvsp) {
  if (!is.na(hgvsp) && hgvsp != "") return(hgvsp)
  annotation <- clean_str(annotation)
  if (is.na(annotation)) return("ne_sais_pas")

  terms <- unlist(strsplit(annotation, "&", fixed = TRUE), use.names = FALSE)

  if ("stop_gained" %in% terms) return("stop_gained")
  if ("frameshift_variant" %in% terms) return("frameshift")
  if ("start_lost" %in% terms) return("start_lost")
  if ("stop_lost" %in% terms) return("stop_lost")
  if ("splice_acceptor_variant" %in% terms || "splice_donor_variant" %in% terms) return("splice_site")
  if ("splice_region_variant" %in% terms) return("splice_region")
  if ("inframe_insertion" %in% terms) return("inframe_insertion")
  if ("inframe_deletion" %in% terms) return("inframe_deletion")
  if ("missense_variant" %in% terms) return("missense")
  if ("synonymous_variant" %in% terms) return("synonymous")
  if ("5_prime_UTR_variant" %in% terms) return("UTR5")
  if ("3_prime_UTR_variant" %in% terms) return("UTR3")
  if ("intron_variant" %in% terms) return("intronique")
  if ("non_coding_transcript_variant" %in% terms ||
      "non_coding_transcript_exon_variant" %in% terms) return("non_coding")
  if ("upstream_gene_variant" %in% terms) return("upstream")
  if ("downstream_gene_variant" %in% terms) return("downstream")
  if ("intergenic_region" %in% terms) return("intergenic")

  "ne_sais_pas"
}

impact_score_fast <- function(impact) {
  impact <- clean_str(impact)
  if (is.na(impact)) return(0L)
  if (impact == "HIGH") return(4L)
  if (impact == "MODERATE") return(3L)
  if (impact == "LOW") return(2L)
  if (impact == "MODIFIER") return(1L)
  0L
}

effect_score_fast <- function(annotation) {
  annotation <- clean_str(annotation)
  if (is.na(annotation)) return(0L)

  terms <- unlist(strsplit(annotation, "&", fixed = TRUE), use.names = FALSE)

  if ("stop_gained" %in% terms) return(100L)
  if ("frameshift_variant" %in% terms) return(95L)
  if ("start_lost" %in% terms) return(90L)
  if ("stop_lost" %in% terms) return(89L)
  if ("splice_acceptor_variant" %in% terms || "splice_donor_variant" %in% terms) return(88L)
  if ("inframe_insertion" %in% terms || "inframe_deletion" %in% terms) return(80L)
  if ("missense_variant" %in% terms) return(70L)
  if ("synonymous_variant" %in% terms || "stop_retained_variant" %in% terms) return(60L)
  if ("splice_region_variant" %in% terms) return(50L)
  if ("5_prime_UTR_variant" %in% terms) return(40L)
  if ("3_prime_UTR_variant" %in% terms) return(39L)
  if ("intron_variant" %in% terms) return(30L)
  if ("non_coding_transcript_variant" %in% terms ||
      "non_coding_transcript_exon_variant" %in% terms) return(20L)
  if ("upstream_gene_variant" %in% terms) return(10L)
  if ("downstream_gene_variant" %in% terms) return(9L)
  if ("intergenic_region" %in% terms) return(5L)

  1L
}

is_coding_biotype_fast <- function(bt) {
  bt <- clean_str(bt)
  if (is.na(bt)) return(FALSE)
  grepl("protein_coding", bt, ignore.case = TRUE)
}

## retourne directement UNE annotation retenue
extract_best_annotation_fast <- function(info_str, gene_str, alt_str) {
  if (is.na(info_str) || is.na(gene_str)) {
    return(c(aa_change_notation = NA_character_,
             variant_annotation_display = "ne_sais_pas"))
  }

  ann_blob <- str_match(info_str, "ANN=([^;]+)")[, 2]
  ann_blob <- clean_str(ann_blob)
  if (is.na(ann_blob)) {
    return(c(aa_change_notation = NA_character_,
             variant_annotation_display = "ne_sais_pas"))
  }

  entries <- strsplit(ann_blob, ",", fixed = TRUE)[[1]]
  if (length(entries) == 0) {
    return(c(aa_change_notation = NA_character_,
             variant_annotation_display = "ne_sais_pas"))
  }

  gene_str <- clean_str(gene_str)
  alt_str  <- clean_str(alt_str)

  best_score <- -Inf
  best_hgvsp <- NA_character_
  best_label <- "ne_sais_pas"

  for (entry in entries) {
    fields <- strsplit(entry, "|", fixed = TRUE)[[1]]

    if (length(fields) < 5) next

    allele       <- clean_str(fields[1])
    annotation   <- if (length(fields) >= 2) clean_str(fields[2]) else NA_character_
    impact       <- if (length(fields) >= 3) clean_str(fields[3]) else NA_character_
    gene_name    <- if (length(fields) >= 4) clean_str(fields[4]) else NA_character_
    gene_id      <- if (length(fields) >= 5) clean_str(fields[5]) else NA_character_
    biotype      <- if (length(fields) >= 8) clean_str(fields[8]) else NA_character_
    hgvsp        <- if (length(fields) >= 11) clean_hgvsp(fields[11]) else NA_character_

    gene_match <- (!is.na(gene_name) && gene_name == gene_str) ||
                  (!is.na(gene_id)   && gene_id   == gene_str)

    if (!gene_match) next

    allele_match <- !is.na(alt_str) && !is.na(allele) && allele == alt_str
    coding_bt    <- is_coding_biotype_fast(biotype)
    has_hgvsp    <- !is.na(hgvsp) && hgvsp != ""

    score <-
      ifelse(gene_match, 1000L, 0L) +
      ifelse(allele_match, 100L, 0L) +
      ifelse(coding_bt, 20L, 0L) +
      ifelse(has_hgvsp, 10L, 0L) +
      effect_score_fast(annotation) +
      impact_score_fast(impact)

    if (score > best_score) {
      best_score <- score
      best_hgvsp <- hgvsp
      best_label <- effect_label_fast(annotation, hgvsp)
    }
  }

  if (!is.finite(best_score)) {
    return(c(aa_change_notation = NA_character_,
             variant_annotation_display = "ne_sais_pas"))
  }

  c(
    aa_change_notation = best_hgvsp,
    variant_annotation_display = best_label
  )
}

## ==== 2. Détection des colonnes d'AF ==== ##
af_cols   <- grep("^[RS]", names(df), value = TRUE)
r_samples <- grep("^R", names(df), value = TRUE)
s_samples <- grep("^S", names(df), value = TRUE)

## ==== 2bis. Colonnes d'impact (compteurs) ==== ##
impact_count_cols <- intersect(c("HIGH", "MODERATE", "MODIFIER", "LOW"), names(df))

## ==== 3. Conversion des colonnes numériques (AF + impacts) ==== ##
df <- df %>%
  mutate(
    across(all_of(af_cols), ~ {
      x <- as.character(.x)
      x[x == "."] <- NA_character_
      x <- str_replace(x, ",", ".")
      as.numeric(x)
    }),
    across(all_of(impact_count_cols), ~ suppressWarnings(as.numeric(.x)))
  )

## ==== 3bis. Ajouter annotation robuste et rapide depuis INFO/ANN ==== ##
if (all(c("INFO", "gene", "variant_alt") %in% names(df))) {
  ann_mat <- t(mapply(
    extract_best_annotation_fast,
    df$INFO,
    df$gene,
    df$variant_alt,
    USE.NAMES = FALSE
  ))

  ann_df <- as.data.frame(ann_mat, stringsAsFactors = FALSE)
  df$aa_change_notation <- ann_df$aa_change_notation
  df$variant_annotation_display <- ann_df$variant_annotation_display
} else {
  df <- df %>%
    mutate(
      aa_change_notation = NA_character_,
      variant_annotation_display = "ne_sais_pas"
    )
}

## ==== 4. Agrégation au niveau gène ==== ##
gene_summary <- df %>%
  group_by(chromosome, gene, gene_start, gene_end) %>%
  summarise(
    ## position milieu du gène
    gene_mid   = (first(gene_start) + first(gene_end)) / 2,

    ## nombre de variants uniques pos/ref/alt
    n_variants = n_distinct(paste(variant_pos, variant_ref, variant_alt, sep = ":")),

    ## longueur du gène (bp)
    gene_len_bp = (first(gene_end) - first(gene_start) + 1),

    ## taux de variants par kb (densité)
    var_per_kb = n_variants / (gene_len_bp / 1000),

    ## AF moyenne par sample sur tous les variants du gène
    across(all_of(af_cols), ~ mean(.x, na.rm = TRUE)),

    ## meilleur impact du gène basé sur les compteurs HIGH/MODERATE/MODIFIER/LOW
    best_impact = dplyr::case_when(
      any(HIGH > 0,     na.rm = TRUE) ~ "HIGH",
      any(MODERATE > 0, na.rm = TRUE) ~ "MODERATE",
      any(MODIFIER > 0, na.rm = TRUE) ~ "MODIFIER",
      any(LOW > 0,      na.rm = TRUE) ~ "LOW",
      TRUE                            ~ "NONE"
    ),

    ## concat de toutes les infos variants pos/ref/alt + annotation choisie
    variant_info = paste(
      unique(paste0(
        variant_pos, " ",
        variant_ref, ">", variant_alt,
        " (", ifelse(is.na(variant_annotation_display) | variant_annotation_display == "",
                     "ne_sais_pas",
                     variant_annotation_display), ")",
        " [H=", coalesce(as.character(HIGH), "NA"),
        ", M=", coalesce(as.character(MODERATE), "NA"),
        ", Md=", coalesce(as.character(MODIFIER), "NA"),
        ", L=", coalesce(as.character(LOW), "NA"), "]"
      )),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    ## moyennes AF R vs S
    mean_R = mean(c_across(all_of(r_samples)), na.rm = TRUE),
    mean_S = mean(c_across(all_of(s_samples)), na.rm = TRUE),
    diff_AF = mean_R - mean_S,

    ## texte AF par sample pour le hover
    af_info = paste(
      paste0(
        af_cols, " = ",
        sprintf("%.3f", c_across(all_of(af_cols)))
      ),
      collapse = "<br>"
    ),

    ## texte complet pour le hover
    hover_text = paste0(
      "Chromosome: ", chromosome,
      "<br>Gene: ", gene,
      "<br>Start: ", gene_start,
      "<br>End: ", gene_end,
      "<br>Mid pos: ", round(gene_mid),
      "<br>Gene length (bp): ", gene_len_bp,
      "<br>#variants: ", n_variants,
      "<br>Variants per kb: ", sprintf("%.3f", var_per_kb),
      "<br>Best impact: ", best_impact,
      "<br>mean_R: ", sprintf("%.3f", mean_R),
      "<br>mean_S: ", sprintf("%.3f", mean_S),
      "<br>diff_AF (R - S): ", sprintf("%.3f", diff_AF),
      "<br><br>AF by sample:<br>", af_info,
      "<br><br>Variants:<br>", variant_info
    )
  ) %>%
  ungroup() %>%
  mutate(
    best_impact = factor(
      best_impact,
      levels = c("HIGH", "MODERATE", "MODIFIER", "LOW", "NONE")
    )
  )

## ==== 4bis. Export table : top 10% des gènes (|diff_AF|) ==== ##
top10_threshold <- quantile(abs(gene_summary$diff_AF), probs = 0.90, na.rm = TRUE)

top10_genes <- gene_summary %>%
  filter(!is.na(diff_AF)) %>%
  mutate(abs_diff_AF = abs(diff_AF)) %>%
  filter(abs_diff_AF >= top10_threshold) %>%
  arrange(desc(abs_diff_AF))

write_tsv(top10_genes, "./7_plotly_gene/top10pct_genes_by_absDiffAF.tsv")

## ==== 5. Génération d'un plotly par chromosome ==== ##
plots_by_chr <- lapply(
  split(gene_summary, gene_summary$chromosome),
  function(dat_chr) {

    layers_order <- c("NONE", "LOW", "MODIFIER", "MODERATE", "HIGH")
    sizeref_chr <- 2 * max(dat_chr$var_per_kb, na.rm = TRUE) / (25^2)

    p <- plot_ly()

    for (imp in layers_order) {
      dat_layer <- dat_chr %>% filter(best_impact == imp)
      if (nrow(dat_layer) == 0) next

      col_rgba <- switch(
        imp,
        "HIGH"     = "rgba(187,11,11,0.95)",
        "MODERATE" = "rgba(255,127,0,0.95)",
        "MODIFIER" = "rgba(199,207,0,0.95)",
        "LOW"      = "rgba(20,148,20,0.95)",
        "NONE"     = "rgba(128,128,128,0.95)",
        "rgba(128,128,128,0.95)"
      )

      p <- p %>%
        add_trace(
          data = dat_layer,
          x = ~gene_mid,
          y = ~diff_AF,
          type = "scatter",
          mode = "markers",
          text = ~hover_text,
          hoverinfo = "text",
          name = imp,
          legendgroup = imp,
          marker = list(
            size = ~var_per_kb,
            sizemode = "area",
            sizeref = sizeref_chr,
            sizemin = 2,
            color = col_rgba,
            line = list(width = 0)
          ),
          showlegend = TRUE
        )
    }

    p <- p %>%
      add_trace(
        data = dat_chr,
        x = ~gene_mid,
        y = ~diff_AF,
        type = "scatter",
        mode = "text",
        text = ~gene,
        textposition = "middle center",
        textfont = list(color = "rgba(0,0,0,0.25)", size = 9),
        hoverinfo = "skip",
        name = "Gene labels",
        showlegend = TRUE,
        visible = "legendonly"
      )

    p %>%
      layout(
        title = unique(dat_chr$chromosome),
        xaxis = list(title = "Gene midpoint (bp)"),
        yaxis = list(title = "mean(AF_R) - mean(AF_S)"),
        legend = list(title = list(text = "Impact / Labels"))
      )
  }
)

## ==== 6. Sauvegarde des plots dans des fichiers HTML ==== ##
if (!rmarkdown::pandoc_available()) {
  stop(
    "Pandoc introuvable : impossible de produire des HTML self-contained.\n",
    "Installe Pandoc (ex: 'sudo apt-get install pandoc' sur Debian/Ubuntu,\n",
    "ou via conda: 'conda install -c conda-forge pandoc'), puis relance."
  )
}

file_names <- character(0)

for (chr in names(plots_by_chr)) {
  p <- plots_by_chr[[chr]]
  safe_chr <- gsub("[^A-Za-z0-9_-]", "_", chr)
  file_name <- paste0("./7_plotly_gene/", safe_chr, ".html")
  saveWidget(p, file = file_name, selfcontained = TRUE)
  file_names <- c(file_names, file_name)
}

cat("Fichiers générés :\n")
cat(paste0(" - ", file_names), sep = "\n")
cat("\nTable exportée :\n")
cat(" - top10pct_genes_by_absDiffAF.tsv\n")
```

## `README.md`

Raw source: lines 1782–1794 of the compilation.

```markdown
# Pipeline d'appel comparatif de variants

## Prérequis
Les fichiers d'appels de variants annotés SnpEff (et leur fichier d'index) ont été générés par échantillon sur une même référence.

## Processus
Placez-les dans 0_samples, dans des dossiers que chacun vous préfixez de la condition, "R_" ou "S_".
Lancez le pipeline avec VarrEffPlot.sh .

## Remarques
Pour personnaliser les graphes, modifier manuellement dans les Rscripts.
Pour relancer une analyse "downstream" (par exemple une figure modifiée), il suffit de supprimer le dossier puis relancer, plutôt que relancer depuis le début.
Utilisez l'option -h (help) pour être guidé et afficher les options notamment de personnalisation des figures statiques ggplot2.
```

## `VarEffPlot.sh`

Raw source: lines 1799–2140 of the compilation.

```bash
#!/bin/bash

args=("$@")

for arg in "${!args[@]}"; do if [[ "${args[arg]}" == "-h" ]]; then
  echo "-> PREREQUISITES"
  echo "   You must have done :"
  echo "   - variant calling with bcftools (mpileup/call) and,"
  echo "   - variant effect prediction annotation with SnpEff and,"
  echo "   - vcf files indexed in .tbi."
  echo "   In folder 0_samples, paste your sample-folders."
  echo "   In each sample-folder, (1) the .vcf.gz and (2) the .vcf.gz.tbi."
  echo "   Folder names must begin by R_ (Resistant sample) or S_ (Susceptible sample)."
  echo "-> DEPENDENCIES"
  echo "   - bcftools installed,"
  echo "   - R with packages ggplot2 and plotly installed."
  echo "-> TO DO : "
  echo "   Then, launch the main script for analysis."
  echo "-> OPTIONS"
  echo "   -l = left limit of genomic position for ggplot (in bp)"
  echo "   -r = right limit of genomic position for ggplot (in bp)"
  echo "   -b = bottom limit of AF difference"
  echo "   -n = add geneID labels on ggplot"
  echo "   -f = force y-axis from 0 to 1 on ggplot"
  exit 0
 fi
done

# 0_tools
mkdir -p "0_tools"

download_tool() {
 file=$1
 url="https://storage.ionos.fr/s/cQNq7E5iQ763Rcw/download?path=%2F&files=$file"
 out="./0_tools/$file"

 if [ ! -s "$out" ]; then
  echo "Downloading $file"
  if command -v curl >/dev/null 2>&1; then
   curl -fL "$url" -o "$out"
  elif command -v wget >/dev/null 2>&1; then
   wget -O "$out" "$url"
  else
   echo "ERREUR: curl ou wget est requis pour télécharger $file" >&2
   exit 1
  fi
 fi
}

download_tool bcftools.tar
download_tool r-plotly.tar
download_tool r-ggplot2.tar

# 1_vcf_reheadered
if [ ! -d "1_vcf_reheadered" ]; then
 mkdir "1_vcf_reheadered"
 for sample in 0_samples/* ; do
  samplename=$(basename "$sample")
  mkdir "./1_vcf_reheadered/$samplename"
  vcf="$(find "$sample" -maxdepth 1 -type f -name "*.vcf.gz" | head -n 1)"
  out="./1_vcf_reheadered/$samplename/$samplename.vcf.gz"
  oldname="$(bcftools query -l "$vcf")"
  mapfile_tmp="$(mktemp)"
  printf "%s\t%s\n" "$oldname" "$samplename" > "$mapfile_tmp"
  bcftools reheader -s "$mapfile_tmp" -o "$out" "$vcf"
  rm "$mapfile_tmp"
  bcftools index -t "$out"
 done
fi

# 2_vcf_merged
if [ ! -d "2_vcf_merged" ]; then
 mkdir "2_vcf_merged"
 bcftools merge -Oz -o ./2_vcf_merged/all.vcf.gz $(find ./1_vcf_reheadered/ -type f -name '*.vcf.gz' -print0 | xargs -0 printf '%s ')
 bcftools index -t ./2_vcf_merged/all.vcf.gz
fi

# 3_gtf_filtered
if [ ! -d "3_gtf_filtered" ]; then
 mkdir "3_gtf_filtered"
 set -euo pipefail
 GTF="$(find ./0_refseq -maxdepth 1 -type f -name '*.gtf' | head -n 1)" # chemin vers ton GTF
 OUTGTF="./3_gtf_filtered/gtf.tsv"

 # TSV: chr, gene_id, type(g/e), start, end
 awk -F'\t' 'BEGIN{OFS="\t"}
   $0 ~ /^#/ {next}
   ($3!="gene" && $3!="exon") {next}
   {
     chrom=$1
     start=$4
     end=$5

     gid="NA"
     if(match($0, /gene_id "[^"]+"/)){
       gid=substr($0, RSTART+9, RLENGTH-10)
     }

     t=($3=="gene" ? "g" : "e")
     print chrom, gid, t, start, end
   }' "$GTF" > "$OUTGTF"

 echo "OK -> $OUTGTF"
fi

# 3_vcf_filtered
if [ ! -d "3_vcf_filtered" ]; then
	mkdir "3_vcf_filtered"
	set -euo pipefail
	# ---- à modifier ----
	IN="./2_vcf_merged/all.vcf.gz"
	OUT="./3_vcf_filtered/vcf.tsv"   # ou 05.final.tsv.gz si tu veux recompresser à la fin
	GTF="$(find ./0_refseq -maxdepth 1 -type f -name '*.gtf' | head -n 1)" # chemin vers ton GTF
	# filtre résistants: AF > 0.4 dans au moins N échantillons R*
	MIN_R_COUNT=1
	# --------------------
	cur="$IN"
	echo "[1] Enlever l'entête '##' (garder #CHROM et variants)"
	step="./3_vcf_filtered/01.no_double_hash.tsv"
	gunzip -c "$cur" | awk 'BEGIN{OFS="\t"} /^##/ {next} {print}' > "$step"
	cur="$step"

	echo "[2] Garder seulement: CHROM POS REF ALT INFO + samples (sans FORMAT)"
	step="./3_vcf_filtered/02.keep_cols.tsv"
	awk 'BEGIN{OFS="\t"}
	     NR==1{
	       # Header VCF: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample...
	       printf "%s\t%s\t%s\t%s\t%s", $1,$2,$4,$5,$8
	       for(i=10;i<=NF;i++) printf "\t%s",$i
	       printf "\n"
	       next
	     }
	     {
	       printf "%s\t%s\t%s\t%s\t%s", $1,$2,$4,$5,$8
	       for(i=10;i<=NF;i++) printf "\t%s",$i
	       printf "\n"
	     }' "$cur" > "$step"
	cur="$step"
	echo "[3] Ajouter HIGH/MODERATE/MODIFIER/LOW (comptés dans INFO) + minDP maxDP (entre samples) + minQUAL maxQUAL=NA + remplacer chaque sample par son AF (AD_alt/DP), manquants -> 0"
	step="./3_vcf_filtered/03.freq_decimal_with_counts_minmax.tsv"
	awk -F'\t' 'BEGIN{OFS="\t"}
	     NR==1{
	       # Colonnes actuelles: #CHROM POS REF ALT INFO sample...
	       printf "%s\t%s\t%s\t%s\t%s\tHIGH\tMODERATE\tMODIFIER\tLOW\tminDP\tmaxDP\tminQUAL\tmaxQUAL", $1,$2,$3,$4,$5
	       for(i=6;i<=NF;i++) printf "\t%s",$i
	       printf "\n"
	       next
	     }
	     {
	       info=$5
	       tmp=info; high     = gsub(/\|HIGH\|/,"",tmp)
	       tmp=info; moderate = gsub(/\|MODERATE\|/,"",tmp)
	       tmp=info; modifier = gsub(/\|MODIFIER\|/,"",tmp)
	       tmp=info; low      = gsub(/\|LOW\|/,"",tmp)

	       mindp=""; maxdp=""
	       for(i=6;i<=NF;i++){
		 s=$i
		 if(s=="./." || s=="." || s==""){
		   $i="0.00"
		   continue
		 }
		 n=split(s,a,":")
		 dp   = (n>=3 ? a[3] : "")
		 ad   = (n>=4 ? a[4] : "")

		 if(dp != "" && dp != "." && dp ~ /^[0-9]+([.][0-9]+)?$/){
		   d = dp + 0
		   if(mindp=="" || d < mindp) mindp=d
		   if(maxdp=="" || d > maxdp) maxdp=d
		 }

		 if(dp=="" || dp=="." || (dp+0)<=0){
		   $i="0.00"
		 } else if(ad=="" || ad=="." ){
		   $i="0.00"
		 } else {
		   split(ad,adp,",")
		   altad = (length(adp)>=2 ? adp[2] : adp[1])
		   if(altad=="" || altad=="." || altad !~ /^[0-9]+([.][0-9]+)?$/){
		     $i="0.00"
		   } else {
		     $i=sprintf("%.2f", (altad+0)/(dp+0))
		   }
		 }
	       }
	       if(mindp=="") mindp="NA"
	       if(maxdp=="") maxdp="NA"
	       minq="NA"; maxq="NA"

	       printf "%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s", $1,$2,$3,$4,$5,high,moderate,modifier,low,mindp,maxdp,minq,maxq
	       for(i=6;i<=NF;i++) printf "\t%s",$i
	       printf "\n"
	     }' "$cur" > "$step"
	cur="$step"
	echo "[4] Supprimer variants si AF > 0.1 dans au moins un S (samples dont le nom commence par 'S')"
	step="./3_vcf_filtered/04.drop_S_gt0.1.tsv"
	awk 'BEGIN{FS=OFS="\t"}
	     NR==1{
	       ns=0
	       for(i=1;i<=NF;i++){
		 if($i ~ /^S/){
		   scol[++ns]=i
		 }
	       }
	       if(ns==0){ print "ERREUR: aucune colonne sample ne commence par S dans l header" > "/dev/stderr"; exit 1 }
	       print; next
	     }
	     {
	       bad=0
	       for(j=1;j<=ns;j++){
		 c=scol[j]
		 if($c!="NA" && ($c+0)>0.1){ bad=1; break }
	       }
	       if(!bad) print
	     }' "$cur" > "$step"
	cur="$step"
	echo "[5] Garder variants si AF > 0.4 dans AU MOINS ${MIN_R_COUNT} R (samples dont le nom commence par 'R')"
	step="./3_vcf_filtered/05.filtered.tsv"
	awk -v MIN_R_COUNT="$MIN_R_COUNT" 'BEGIN{FS=OFS="\t"}
	     NR==1{
	       nr=0
	       for(i=1;i<=NF;i++){
		 if($i ~ /^R/){
		   rcol[++nr]=i
		 }
	       }
	       if(nr==0){ print "ERREUR: aucune colonne sample ne commence par R dans l header" > "/dev/stderr"; exit 1 }
	       print; next
	     }
	     {
	       pass=0
	       for(j=1;j<=nr;j++){
		 c=rcol[j]
		 if($c!="NA" && ($c+0)>0.4){ pass++ }
		 if(pass>=MIN_R_COUNT){ break }
	       }
	       if(pass>=MIN_R_COUNT) print
	     }' "$cur" > "$step"
	cur="$step"
	echo "[6] Ajouter les HCON_IDs (gènes contenant le variant) en 2ème colonne (après CHR)"
	step="./3_vcf_filtered/06.add_gene_from_gtf.tsv"
	awk -F'\t' 'BEGIN{OFS="\t"}
	# --- 1) lecture du GTF (premier fichier) ---
	FNR==NR{
	  if($0 ~ /^#/) next
	  if($3!="gene") next
	  chrom=$1
	  start=$4+0
	  end=$5+0
	  gid=""
	  if(match($0, /gene_id "[^"]+"/)){
	    gid=substr($0, RSTART+9, RLENGTH-10)
	  }
	  if(gid=="") next
	  k = chrom SUBSEP (++n[chrom])
	  s[k]=start
	  e[k]=end
	  id[k]=gid
	  next
	}
	# --- 2) lecture du TSV ---
	FNR==1{
	  # header attendu: CHROM POS REF ALT INFO ...
	  # sortie: CHROM HCON_ID POS REF ALT INFO ...
	  printf "%s\tHCON_ID\t", $1
	  for(i=2;i<=NF;i++){
	    printf "%s", $i
	    if(i<NF) printf "\t"
	  }
	  printf "\n"
	  next
	}
	{
	  chrom=$1
	  pos=$2+0

	  out="NA"
	  if(chrom in n){
	    j = ptr[chrom]; if(j=="") j=1
	    while(j <= n[chrom] && e[chrom SUBSEP j] < pos) j++
	    ptr[chrom]=j
	    tmp=""
	    k=j
	    while(k <= n[chrom] && s[chrom SUBSEP k] <= pos){
	      if(pos <= e[chrom SUBSEP k]){
	        gid=id[chrom SUBSEP k]
	        if(tmp=="") tmp=gid
	        else tmp=tmp","gid
	      }
	      k++
	    }
	    if(tmp!="") out=tmp
	  }

	  # print reorder: CHROM, HCON_ID, POS, REF, ALT, INFO, HIGH..., samples...
	  printf "%s\t%s\t", $1, out
	  for(i=2;i<=NF;i++){
	    printf "%s", $i
	    if(i<NF) printf "\t"
	  }
	  printf "\n"
	}' "$GTF" "$cur" > "$OUT"
	echo "OK -> $OUT"
fi

if [ ! -d "4_join" ]; then
 mkdir "4_join"
 Rscript ./0_Rscripts/4_gtf_vcf_join.r
fi

if [ ! -d "4_gene_var_eff" ]; then
 mkdir "4_gene_var_eff"
 Rscript ./0_Rscripts/4_variants_to_gene_effect_table.r
fi

if [ ! -d "5_var_density" ]; then
 mkdir "5_var_density"
 Rscript ./0_Rscripts/5_var_density.r
fi

if [ ! -d "6_plotly_var" ]; then
 mkdir "6_plotly_var"
 Rscript ./0_Rscripts/6_plotly_var.r
fi

if [ ! -d "7_plotly_gene" ]; then
 mkdir "7_plotly_gene"
 Rscript ./0_Rscripts/7_plotly_gene.r
fi

if [ ! -d "7_ggplot_gene" ]; then
 mkdir "7_ggplot_gene"
 for arg in "${!args[@]}"; do if [[ "${args[arg]}" == "-l" ]]; then l=${args[arg+1]}; fi; done
 for arg in "${!args[@]}"; do if [[ "${args[arg]}" == "-r" ]]; then r=${args[arg+1]}; fi; done
 for arg in "${!args[@]}"; do if [[ "${args[arg]}" == "-b" ]]; then b=${args[arg+1]}; fi; done
 n=""
 for arg in "${!args[@]}"; do if [[ "${args[arg]}" == "-n" ]]; then n="-n"; fi; done
 f=""
 for arg in "${!args[@]}"; do if [[ "${args[arg]}" == "-f" ]]; then f="-f"; fi; done
 Rscript ./0_Rscripts/7_ggplot_gene.r -l ${l} -r ${r} -b ${b} ${n} ${f}
fi
```
