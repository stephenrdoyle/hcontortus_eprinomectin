# Documentation du pipeline transcriptomique et de l’analyse différentielle

> **Objet du document.** Cette page transforme un cahier de laboratoire brut en une documentation scientifique destinée à accompagner le dépôt du code lié à une publication. Elle décrit les **étapes, commandes, paramètres, outils, versions, fichiers, conventions et états d’activation consignés dans l’instantané fourni**.  
>
> L’objectif n’est pas de présenter le dépôt comme un workflow immédiatement réexécutable, mais de permettre au lecteur de comprendre précisément **ce qui a été codé et comment les analyses ont été conduites ou préparées**.

## Métadonnées générales

| Champ | Information consignée |
|---|---|
| Organisme | *Haemonchus contortus* |
| Données | RNA-seq Illumina paired-end |
| Analyse principale | alignement génomique, comptage génique et expression différentielle |
| Génome de référence | `haemonchus_contortus.PRJEB506.WBPS18.genomic.fa` |
| Annotation principale | `haemonchus_contortus.PRJEB506.WBPS18.canonical_geneset.gtf` |
| Transcriptome alternatif | `haemonchus_contortus.PRJEB506.WBPS18.mRNA_transcripts.fa` |
| Infrastructure principale | cluster HPC avec environnement Slurm et modules logiciels |
| Infrastructure alternative | images Singularity pour FastQC, MultiQC, Cutadapt, HISAT2 et samtools |
| Version de R appelée | R 4.3.2 |
| Auteur indiqué dans le README historique | Robin Lioutaud |
| Unité indiquée | UMR INRAE/ENVT 1436 InTheRes, Toulouse, France |
| Licence indiquée | CeCILL |
| Publication associée | À renseigner |
| DOI / accession des lectures | À renseigner |
| Date exacte des analyses | Non consignée |
| Commit correspondant aux résultats publiés | À renseigner |

## Portée de cette compilation

Le cahier contient un pipeline RNA-seq allant du contrôle qualité jusqu’aux visualisations de l’expression différentielle. Il comprend :

- un chemin principal utilisant **FastQC → MultiQC → Cutadapt → HISAT2/samtools → Sambamba → featureCounts → DESeq2 → volcano plots** ;
- une branche alternative de pseudo-alignement **Kallisto → Sleuth** ;
- une analyse alternative **edgeR** ;
- plusieurs normalisations ou transformations exploratoires : **CPM, MRM/DESeq2, TMM, TPM et RPKM** ;
- une exploration des échantillons par **PCA** et heatmap de distances ;
- des essais archivés avec **HTSeq** et **MMQuant**.

Aucun script d’enrichissement fonctionnel, de GSEA, de variant calling ou de fusion avec une base d’annotation n’est présent dans cette compilation, même si le README historique mentionne un « GSEA manuel ».

### Interprétation des blocs conditionnels

Les scripts ont été développés de manière itérative :

- `if (TRUE)` ou `if true` : bloc actif dans l’instantané ;
- `if (FALSE)` ou `if false` : bloc désactivé dans l’instantané ;
- un bloc désactivé peut correspondre à une étape exécutée antérieurement lorsque ses sorties sont relues ensuite ;
- l’état du code ne constitue donc pas, à lui seul, une preuve d’exécution pour chaque résultat de la publication.

## Vue d’ensemble

```text
FASTQ paired-end bruts
        │
        ├── FastQC ── MultiQC
        │
        └── Cutadapt
                │
                └── HISAT2 + samtools
                        │
                        └── BAM triés et indexés
                                │
                                └── Sambamba markdup
                                        │
                                        └── featureCounts
                                                │
                                                ├── CPM / MRM / TMM / TPM
                                                ├── PCA + heatmap de distances
                                                └── DESeq2
                                                        │
                                                        └── volcano plots Plotly

Branche alternative :
FASTQ → Kallisto → abondances transcriptomiques → Sleuth
```

## Arborescence analytique documentée

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

Le répertoire `.git/` figurait également dans la compilation brute. Ses fichiers internes et les hooks d’exemple fournis par Git ne sont pas des étapes analytiques et ne sont pas reproduits dans l’annexe.

## État des étapes dans l’instantané

| Répertoire | Fonction | État visible |
|---|---|---|
| `01_fastqc` | contrôle qualité initial | actif |
| `02_multiqc` | agrégation des rapports | actif |
| `03_cutadapt` | retrait des adaptateurs | actif |
| `04_hisat2` | index, alignement, tri, index BAM et métriques | actifs |
| `04_hisat2/autre/kallisto` | quantification alternative | index et quantification désactivés ; concaténation des TSV active |
| `05_MarkDuplicates` | marquage des duplicats | actif |
| `05_MarkDuplicates/quality.sh` | métriques post-duplicats | génération `flagstat/stats` désactivée ; synthèse CSV active |
| `11_featurecounts` | comptage génique | calcul featureCounts désactivé ; renommage d’un `cnt.csv` existant actif |
| `13_normalisation_*` | CPM, MRM, TMM et TPM | actif lorsque chaque script est lancé |
| `14_normalisation` | brouillon RPKM et choix de normalisation | calcul RPKM actif dans le brouillon ; lanceur configuré pour `tpm.r`, fichier absent de cette compilation |
| `15_exploration` | PCA et heatmap de distances | fonction active ; sélection du fichier d’entrée à vérifier |
| `16_deseq2` | analyse différentielle principale | active |
| `16_deseq2/autres/edgeR` | analyse alternative multifactorielle | script autonome |
| `16_deseq2/autres/sleuth` | exploration des résultats Kallisto | script autonome |
| `17_DGEplot` | volcano plots | actif ; heatmap DGE désactivée |
| `main.sh` | orchestration générale | préparée mais incompatible avec les noms actuels des dossiers |

# Organisation des données

Le README historique prévoit deux répertoires frères :

```text
<racine>/
├── RAW/   # FASTQ, adaptateurs, génome et annotations
└── RES/   # dépôt du pipeline et résultats par étape
```

## Entrées explicitement référencées

| Catégorie | Noms ou motifs utilisés |
|---|---|
| FASTQ bruts | `../../RAW/*.fastq*` |
| Chemin absolu des FASTQ | `/work/user/rlioutaud/RAW/hcon/rnaseq` |
| Paires de lectures | `<sample>_1.fastq` et `<sample>_2.fastq` |
| Adaptateurs | `../../RAW/adapter1*.txt`, `../../RAW/adapter2*.txt` |
| Génome | `*.genomic.fa*` |
| Index FASTA | `haemonchus_contortus.PRJEB506.WBPS18.genomic.fa.fai` |
| Annotation GTF | `haemonchus_contortus.PRJEB506.WBPS18.canonical_geneset.gtf` |
| Annotation GFF3 | `haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3` |
| Transcriptome Kallisto | `haemonchus_contortus.PRJEB506.WBPS18.mRNA_transcripts.fa` |
| Association lane–échantillon | `lane_sample_IDs.list` |
| Métadonnées Kallisto/Sleuth | `rnaseq_metadata.txt` |
| Correspondance annotation–fonction | `AnnotationToGeneFunction.csv` |
| Comptages principaux | `11_featurecounts/cnt.csv` |

## Convention de nommage des échantillons

Le code exploite des noms séparés par des points, par exemple :

```text
T.R.ARA.F.1
```

| Position ou code | Valeurs observées | Interprétation utilisée |
|---|---|---|
| `R` / `S` | résistance | Resistant / Susceptible |
| `ARA`, `BET`, `BUN`, `MOU`, `CHI`, `LUC` | population ou souche | groupes biologiques |
| `F` / `M` | sexe | Female / Male |
| entier terminal | réplicat | numéro de réplicat |
| premier élément, par exemple `T` | non défini | signification non explicitée |

Le README donne l’exemple `R.LOC.M.1.fq.gz`, mais les scripts de prétraitement utilisent principalement les suffixes `_1.fastq` et `_2.fastq`.

# Logiciels, environnements et versions

| Logiciel / environnement | Version consignée | Mode d’appel | Remarque |
|---|---:|---|---|
| Ubuntu | 20.04 | base des définitions Singularity | conteneurs FastQC, MultiQC, Cutadapt et HISAT2 |
| FastQC | 0.12.1 | module HPC | conteneur : 0.11.9 |
| MultiQC | 1.14 | module HPC | installation `pip3 install multiqc` non figée dans le conteneur |
| Cutadapt | 4.3 | module HPC | installation `pip3 install cutadapt` non figée dans le conteneur |
| HISAT2 | 2.2.1 | module HPC | conteneur : 2.2.0 |
| samtools | 1.19 | module HPC | utilisé pour conversion, tri, index et métriques |
| Python | 3.11.1 | module HPC | chargé dans l’étape HISAT2 |
| Kallisto | 0.50.1 | module HPC | branche alternative |
| Sambamba | 1.0.1 | module HPC | marquage des duplicats |
| MMQuant | 1.0.9 | module HPC | méthode de comptage alternative |
| R | 4.3.2 | exécutable absolu | chemins `/save/user/...` et `/work/user/...` |
| Singularity | non consignée | images `.sif` | version du moteur absente du cahier |
| Rsubread / featureCounts | non consignée | package R | version du package absente |
| DESeq2 | non consignée | package R | analyse principale et normalisation MRM |
| edgeR | non consignée | package R | TMM et analyse alternative |
| Sleuth | non consignée | package R | branche Kallisto |
| ggplot2 | non consignée | package R | graphiques |
| Plotly | non consignée | package R | PCA et volcano plots interactifs |
| pheatmap | non consignée | package R | heatmaps |
| tidyverse | non consignée | package R | manipulation de données |
| RColorBrewer | non consignée | package R | palettes de heatmap |
| patchwork | non consignée | package R | assemblage des PCA Sleuth |

Aucun `sessionInfo()`, fichier `renv.lock`, environnement Conda ou gel des dépendances Python n’est présent.

# Pipeline principal

## 1. Contrôle qualité initial — FastQC

**Répertoire :** `01_fastqc/`  
**Entrées :** FASTQ paired-end dans `../../RAW/`  
**Sorties :** rapports FastQC dans le répertoire courant.

### HPC

```bash
module load bioinfo/FastQC/0.12.1

fastqc -t 4 -o ./ ../../RAW/${file}_1.fastq* &
fastqc -t 4 -o ./ ../../RAW/${file}_2.fastq*
```

Les identifiants sont construits en supprimant tout ce qui suit le premier underscore :

```bash
FILENAME=$(echo "$FILENAME" | sed 's/_.*//')
```

Les échantillons sont divisés en groupes de six. Pour chaque échantillon, les deux mates sont traités dans un sous-processus lancé en arrière-plan. Le maximum théorique est de douze processus FastQC à quatre threads, soit 48 threads.

### Singularity

```bash
singularity exec ../fastqc.sif fastqc -t 4 -o ./ ../../RAW/${file}_1.fastq*
singularity exec ../fastqc.sif fastqc -t 4 -o ./ ../../RAW/${file}_2.fastq*
```

La définition du conteneur installe FastQC 0.11.9 sur Ubuntu 20.04, avec Java 11.

## 2. Agrégation des contrôles qualité — MultiQC

**Répertoire :** `02_multiqc/`

```bash
module load bioinfo/MultiQC/1.14
multiqc --outdir ./ ../*fastqc
```

Alternative :

```bash
singularity exec ../multiqc.sif multiqc --outdir ./ ../*fastqc
```

La version installée par `pip3 install multiqc` dans la définition Singularity n’est pas fixée.

## 3. Retrait des adaptateurs — Cutadapt

**Répertoire :** `03_cutadapt/`  
**Entrées :** FASTQ paired-end et fichiers d’adaptateurs  
**Sorties :** `<sample>_1.trim.fastq`, `<sample>_2.trim.fastq`.

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

Aucun seuil explicite de qualité, de longueur minimale, de taux d’erreur ou de nombre de bases ambiguës n’est fourni. Sur HPC, toutes les paires sont lancées en arrière-plan sans limite de concurrence explicite, malgré `#SBATCH --ntasks=1`.

## 4. Alignement — HISAT2 et samtools

**Répertoire :** `04_hisat2/`  
**Entrées :** génome et FASTQ trimés  
**Sorties :** SAM temporaires, BAM triés, index BAM et métriques.

### Détection des chemins

```bash
PathGenome=$(find "../../RAW/" -name "*.genomic.fa*" -print -quit)
Path=$(find ../../RES -type d -name "*cutadapt*" -print)
```

Si plusieurs dossiers correspondent à `*cutadapt*`, la variable `Path` peut contenir plusieurs lignes.

### Construction de l’index

```bash
IndexName="hisat2idx"
hisat2-build -p 8 $PathGenome $IndexName
chmod 777 hisat2idx*
```

### Alignement paired-end et read groups

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

Aucune option spécifique de strandedness ou de transcriptome connu (`--rna-strandness`, fichier de splice sites ou exons) n’est consignée.

### Conversion, tri, index et métriques

```bash
samtools view -bS "${fastq}.sam" |
  samtools sort -@ 4 -o "${fastq}.sorted.bam"

rm "${fastq}.sam"
samtools index -@ 4 "$bam"
samtools flagstat "$BAM_FILE" > "${BASENAME}_flagstat.txt"
samtools stats "$BAM_FILE" > "${BASENAME}_stats.txt"
```

Quatre échantillons sont lancés simultanément et chaque HISAT2 reçoit huit threads. Le besoin théorique atteint donc 32 threads, alors que le commentaire Slurm indique `--ntasks=8`.

## 5. Marquage des duplicats — Sambamba

**Répertoire :** `05_MarkDuplicates/`  
**Entrées :** `*.sorted.bam`  
**Sorties :** `*.mkdup.sorted.bam`.

```bash
module load bioinfo/Sambamba/1.0.1

sambamba markdup -t=8 \
  "${input}/${bam}.sorted.bam" \
  "./${bam}.mkdup.sorted.bam"
```

Trois BAM sont traités en parallèle, soit jusqu’à 24 threads. Les duplicats sont **marqués**, pas supprimés. Aucune indexation des BAM post-marquage n’est explicitement consignée.

### Synthèse de qualité

Les blocs suivants sont désactivés :

```bash
samtools flagstat "$BAM_FILE" > "${BASENAME}_flagstat.txt"
samtools stats "$BAM_FILE" > "${BASENAME}_stats.txt"
```

Le bloc actif lit pourtant ces fichiers et écrit `summary_metrics.csv` :

```text
File,Total Reads,Mapped Reads,Alignment Rate,
Uniquely Mapped Reads,Duplication Rate,Coverage Depth
```

Calculs principaux :

```bash
alignment_rate = mapped_reads / total_reads * 100
uniquely_mapped_reads = (mapped_reads - secondary) / total_reads * 100
duplication_rate = duplicates / total_reads * 100
coverage_depth = bases_mapped_cigar / genome_size
```

Le terme `Uniquely Mapped Reads` est une approximation interne fondée sur la soustraction des alignements secondaires, et non la métrique « uniquely mapped » produite par un aligneur. Les variables nommées `rrna_rate` et `exonic_rate` sont calculées à partir de lignes `QC-passed reads` et `QC-failed reads`, mais ne sont pas exportées et ne mesurent pas directement le rRNA ou l’exonicité.

## 6. Comptage génique — featureCounts

**Répertoire :** `11_featurecounts/`  
**Entrée attendue :** BAM post-marquage  
**Sortie :** `cnt.csv`.

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

Export :

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

Le calcul est désactivé (`if (FALSE)`). Le bloc actif relit un `cnt.csv` existant et raccourcit les noms d’échantillons en supprimant le préfixe et le suffixe communs.

Le motif `.*MarkDuplicates.*\\.bam$` doit être vérifié : si `MarkDuplicates` apparaît uniquement dans le nom du dossier et non dans le nom du BAM, la recherche peut ne sélectionner aucun fichier.

Aucun paramètre explicite de strandedness, de type de feature ou d’attribut GTF n’est fourni ; les valeurs par défaut de `featureCounts` s’appliquent aux options absentes.

### Alternatives archivées

#### HTSeq

Commande modèle entièrement commentée :

```bash
htseq-count -f bam -r pos -s no -i gene_id output.bam reference.gtf > counts.txt
```

#### MMQuant 1.0.9

```bash
module load bioinfo/mmquant/1.0.9
mmquant --bam "$bam" --annotation "$annotation" --output "${filename}_counts.csv"
```

Les chemins de cette branche (`../../2xAlignement/hisat2/`, `../featureCounts/annotations.gtf`) appartiennent à une ancienne organisation.

## 7. Quantification alternative — Kallisto

**Répertoire :** `04_hisat2/autre/kallisto/`

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

L’index et la quantification sont désactivés dans l’instantané. La phase active extrait la quatrième colonne des `abundance.tsv`, correspondant normalement à `est_counts`, puis colle les colonnes dans `counts.tsv` :

```bash
cut -f4 -d$'\t' "$tsv" | tail -n +2 > "${filename}.part.tsv"
paste *.part.tsv > counts.tsv
```

Les identifiants de transcrits ne sont pas conservés dans cette concaténation. L’ordre identique des lignes entre tous les fichiers est donc implicitement supposé.

## 8. Normalisations et tables exploratoires

Toutes les méthodes partent d’un `cnt.csv` comprenant `GeneID`, `Length` et les comptages.

| Dossier | Méthode | Sortie | Formule ou fonction |
|---|---|---|---|
| `13_normalisation_cpm` | CPM | `cnt.cpm.csv` | `count / total_mapped × 10^6` |
| `13_normalisation_mrm` | facteurs de taille DESeq2 | `cnt.mrm.csv` | `estimateSizeFactors()` puis `counts(normalized=TRUE)` |
| `13_normalisation_tmm` | TMM edgeR « déCPMisé » | `cnt.tmm.csv` | TMM-CPM puis multiplication par `TotalMapped / 10^6` |
| `13_normalisation_tpm` | TPM manuel | `cnt.tpm.csv` | division par longueur puis somme ramenée à `10^6` |
| `14_normalisation` | RPKM, brouillon | `rpkm.csv` ou `tmm_rpkm.csv` | CPM × `1000 / Length` |

La normalisation MRM crée un objet DESeq2 avec un design `~1`, uniquement pour estimer les facteurs de taille. L’analyse différentielle principale repart des **comptages bruts** et estime elle-même ses facteurs.

Le script `14_normalisation/sh.sh` appelle `tpm.r`, alors que ce fichier n’apparaît pas dans la compilation. Les appels à `tmm.r`, `cpm.r` et `mrm.r` sont commentés.

## 9. Exploration des échantillons — PCA et heatmap

**Répertoire :** `15_exploration/`  
**Sorties visées :** `plots.html`, `ClusteredHeatmap.svg`.

Étapes :

1. lecture d’une table normalisée ;
2. retrait de la colonne de longueur ;
3. transposition échantillons × gènes ;
4. ajout de `epsilon = 1` ;
5. transformation optionnelle ;
6. suppression des gènes constants ;
7. `prcomp(..., scale. = TRUE)` ;
8. PCA Plotly colorées selon chaque composante du nom d’échantillon ;
9. distance `1 - cor` et heatmap hiérarchique.

Transformations codées :

```r
if (tf == "rlog") {
  df <- log(df)
} else if (tf == "vst") {
  variances <- apply(df, 2, var)
  vst_transform <- sqrt(variances)
  df <- log(df + 1) / vst_transform
}
```

Ces transformations ne sont **pas** les fonctions `rlog()` et `vst()` de DESeq2.

Le pourcentage de variance est recalculé ainsi :

```r
varexp <- summary(pca)$importance[2, ]
varexp <- varexp^2
varexp <- varexp / sum(varexp) * 100
```

La proportion de variance fournie par `prcomp` est donc élevée au carré une seconde fois ; les valeurs affichées ne correspondent pas aux proportions conventionnelles.

L’appel final construit un chemin littéral :

```r
file.path(
  list.files("../", pattern = ".*normalisation_mrm.*", full.names = TRUE),
  ".*mrm.csv"
)
```

Il ne recherche pas réellement les fichiers avec une expression régulière à l’intérieur du dossier. Ce point doit être vérifié par rapport à l’exécution historique.

## 10. Analyse différentielle principale — DESeq2

**Répertoire :** `16_deseq2/`  
**Entrée :** `../11_featurecounts/cnt.csv`  
**Sortie principale :** `cnt.deseq2.dge.csv`.

### Préparation

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

### Contraste sélectionné

Le code actif conserve les femelles :

```r
counts <- SeleCols(counts, "\\.F\\.")
```

Référence :

```r
reference <- SeleCols(reference, "\\.CHI\\.|\\.LUC\\.")
```

Comparaison :

```r
comparison <- SeleCols(comparison, "\\.BET\\.")
```

Le contraste est donc :

```text
BET femelles versus CHI/LUC femelles
```

sous réserve que les noms de colonnes suivent exactement les conventions attendues.

### Modèle et préfiltrage

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

### Normalisation et test

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

Avec les niveaux `ref`, puis `comp`, le coefficient par défaut représente normalement `comp` relativement à `ref`. Un `log2FoldChange` positif correspond donc à une expression plus élevée dans BET femelles que dans CHI/LUC femelles.

Le script exporte temporairement les comptages normalisés vers `11_featurecounts/cnt.mrm.csv`, puis supprime tous les fichiers du dossier source dont le nom correspond à `.mrm.csv`.

### Figures de dispersion

```r
pdf("cnt_dispersion.pdf")
plotDispEsts(dds)
dev.off()
```

Un second bloc tente de produire un SVG personnalisé à partir de `dispersionFunction(dds)`. Cet objet est une fonction de dispersion, et ce graphique ne doit pas être assimilé sans vérification au tracé standard de `plotDispEsts()`.

Le MA plot est présent mais désactivé.

## 11. Analyse alternative — edgeR

**Répertoire :** `16_deseq2/autres/edgeR/`

Le script lit une ancienne table :

```r
/save/user/rlioutaud/pipelines/rnaseq/3vCounting/featureCounts/counts.txt
```

Les positions 2 à 4 des noms d’échantillons sont converties en facteurs, puis le modèle suivant est construit :

```r
dge <- DGEList(counts = counts, group = coldata$RS)
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + RS + FM + RS:FM, data = coldata)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
results <- topTags(qlf)
```

Points à noter :

- le commentaire initial mentionne DESeq2, mais le package utilisé est edgeR ;
- aucune étape explicite `estimateDisp()` n’est présente avant `glmQLFit()` ;
- `coef=2` désigne la deuxième colonne de la matrice de design, sans nommer explicitement le contraste ;
- seuls les résultats affichés dans `log.log` sont consignés ; aucun export CSV n’est effectué.

Cette branche est à considérer comme exploratoire.

## 12. Branche Kallisto/Sleuth

**Répertoire :** `16_deseq2/autres/sleuth/`

```r
metadata <- read.table("rnaseq_metadata.txt", header = TRUE)
so <- sleuth_prep(
  metadata,
  extra_bootstrap_summary = TRUE,
  num_cores = 2
)
```

Le script :

- extrait une matrice d’`est_counts` ;
- réalise une PCA ;
- joint les scores aux covariables ;
- génère trois PCA colorées par `population`, `sex` et `drug` ;
- calcule une matrice de divergence de Jensen–Shannon à partir des TPM ;
- produit une heatmap annotée.

Sorties :

```text
figure_rnaseq_pca_pop_sex_drug.png
figure_rnaseq_heatmap_pop_sex_drug.pdf
figure_rnaseq_heatmap_pop_sex_drug.png
```

Aucun modèle différentiel Sleuth (`sleuth_fit`, `sleuth_wt` ou `sleuth_lrt`) n’est consigné. Cette branche réalise donc surtout une exploration des échantillons.

## 13. Volcano plots — Plotly

**Répertoire :** `17_DGEplot/`  
**Entrées :** `../16_deseq2/*.dge.csv`  
**Sorties :** `*_volcano_plot.html`.

```r
names(dge)[names(dge) == "log2FoldChange"] <- "log2FC"
names(dge)[names(dge) == "padj"] <- "neglog10padj"
dge$neglog10padj <- -log10(dge$neglog10padj)
```

Seuils :

| Classe | Critère |
|---|---|
| non significatif ou effet faible | `abs(log2FC) < 1` ou `-log10(padj) < 1.3` |
| sous-exprimé | `log2FC <= -1` et `-log10(padj) >= 1.3` |
| sur-exprimé | `log2FC >= 1` et `-log10(padj) >= 1.3` |

`1.3` correspond approximativement à `padj = 0,05`.

```r
htmlwidgets::saveWidget(
  volcan,
  paste0(filename, "_volcano_plot.html")
)
```

La heatmap DGE est désactivée. Elle dépend en outre d’objets `dds` et `res` non créés dans ce script et de fichiers externes `Mcounts_normalized.csv` et `Fcounts_normalized.csv`.

# Inventaire des sorties

| Étape | Sorties principales |
|---|---|
| FastQC | rapports HTML et archives ZIP |
| MultiQC | rapport MultiQC |
| Cutadapt | `*_1.trim.fastq`, `*_2.trim.fastq` |
| HISAT2/samtools | `*.sorted.bam`, `.bai`, `*_flagstat.txt`, `*_stats.txt` |
| Sambamba | `*.mkdup.sorted.bam` |
| Qualité post-duplicats | `summary_metrics.csv` |
| featureCounts | `cnt.csv` |
| CPM | `cnt.cpm.csv` |
| MRM | `cnt.mrm.csv` |
| TMM | `cnt.tmm.csv` |
| TPM | `cnt.tpm.csv` |
| RPKM exploratoire | `rpkm.csv`, `tmm_rpkm.csv` |
| PCA/heatmap | `plots.html`, `ClusteredHeatmap.svg` |
| DESeq2 | `cnt.deseq2.dge.csv`, `cnt_dispersion.pdf`, `cnt_dispersion.svg`, `log.log` |
| Volcano plot | `cnt.deseq2.dge_volcano_plot.html` |
| Kallisto | `abundance.tsv`, `abundance.h5`, `run_info.json`, `counts.tsv` |
| Sleuth | PCA PNG, heatmap PDF et PNG |

# Points de vigilance

1. **Versions divergentes entre HPC et conteneurs.** FastQC 0.12.1/0.11.9 et HISAT2 2.2.1/2.2.0.
2. **Dépendances non figées.** MultiQC et Cutadapt sont installés par `pip` sans version ; les packages R n’ont pas de versions consignées.
3. **Chemins absolus.** Plusieurs scripts sont liés aux comptes `/work/user/rlioutaud` et `/save/user/rlioutaud`.
4. **Directives Slurm mal positionnées.** Les lignes `#SBATCH` apparaissent après une instruction shell ; elles ne sont normalement interprétées par Slurm que lorsqu’elles précèdent toute instruction exécutable.
5. **Parallélisation non alignée sur les ressources.** HISAT2 peut utiliser 32 threads pour une demande commentée de huit tâches ; Cutadapt lance tous les échantillons simultanément.
6. **Permissions très ouvertes.** `chmod 777` est appliqué à de nombreux fichiers.
7. **Identifiants tronqués au premier underscore.** `sed 's/_.*//'` peut fusionner ou tronquer des identifiants complexes.
8. **featureCounts n’est pas recalculé dans l’instantané.** Le script actif suppose `cnt.csv` déjà présent.
9. **Sélection des BAM à vérifier.** Le motif featureCounts inclut `MarkDuplicates` dans le motif de fichier.
10. **Strandedness absente.** Aucune orientation de librairie n’est consignée pour HISAT2 ou featureCounts.
11. **QC post-duplicats incohérent dans l’état courant.** La synthèse active dépend de métriques dont la génération est désactivée.
12. **Concaténation Kallisto sans identifiants.** `counts.tsv` repose uniquement sur l’ordre des lignes.
13. **Transformations dites rlog/vst non standard.** Elles ne correspondent pas aux transformations DESeq2.
14. **Pourcentage de variance PCA non standard.** Les proportions sont élevées au carré deux fois.
15. **Contraste DESeq2 codé en dur.** Seules les femelles BET versus CHI/LUC sont testées.
16. **Sortie edgeR incomplète.** Pas d’estimation de dispersion explicite ni de CSV final.
17. **Heatmap DGE non autonome.** Le bloc désactivé dépend d’objets et fichiers externes.
18. **Orchestration générale non fonctionnelle telle quelle.** `main.sh` ne cible pas les dossiers `01_fastqc`, `02_multiqc`, etc.
19. **GSEA annoncé mais absent.** Aucun script de GSEA n’est présent dans cette compilation.
20. **Date et provenance exacte des résultats non consignées.** L’instantané doit être relié manuellement au commit et aux figures de la publication.

# Exécution globale et dépôt Git

Le script `main.sh` recherche les dossiers avec :

```bash
if [[ $folder =~ ^[0-9]*v ]]; then
```

Les dossiers réels sont nommés `01_fastqc`, `02_multiqc`, etc. Ils ne contiennent pas la lettre `v` après leur préfixe numérique. Le script ne les sélectionne donc pas dans son état actuel.

Le README historique indique explicitement que le pipeline n’était pas finalisé et que les étapes devaient être exécutées manuellement. Il mentionne également un risque de retours chariot Windows :

```bash
sed -i 's/\r//' LeScriptProblematique.sh
```

# Informations à compléter avant publication

| Élément | Statut |
|---|---|
| DOI de l’article | à renseigner |
| accession SRA/ENA des lectures | à renseigner |
| tableau des échantillons | à ajouter |
| date ou période d’exécution | à ajouter si disponible |
| commit Git associé aux résultats | à fixer |
| contraste exact présenté dans l’article | à confirmer |
| fichiers ayant alimenté chaque figure | à lister |
| rapports FastQC/MultiQC | à archiver ou lier |
| versions des packages R | à fournir uniquement si historiquement récupérables |
| sommes de contrôle des références | à ajouter si disponibles |

Il est préférable de conserver les scripts historiques tels quels dans un sous-répertoire `scripts/` et d’utiliser ce document comme description interprétée. Les incohérences ne doivent pas être corrigées silencieusement dans la documentation d’une analyse déjà réalisée.

# Annexe — fichiers analytiques reproduits sans modification

Cette annexe reproduit les fichiers analytiques, les lanceurs, les définitions de conteneurs, `main.sh` et le README historique. Les fichiers internes de `.git/` sont exclus.


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
