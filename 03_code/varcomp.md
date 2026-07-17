# Documentation du pipeline comparatif de variants

> **Objet du document.** Cette page transforme un cahier de laboratoire brut en une documentation scientifique destinée à accompagner un dépôt Git associé à une publication. Elle décrit les **commandes, paramètres, outils, fichiers, filtres, calculs et visualisations effectivement consignés dans le code fourni**.
>
> Le but n’est pas de présenter le dépôt comme un workflow entièrement reproductible ou immédiatement réanalysable. Le but est de permettre à un lecteur de savoir précisément **quelles opérations ont été programmées et selon quelles règles**.

## Résumé

Ce dépôt met en œuvre une analyse comparative de variants entre deux groupes d’échantillons :

- échantillons **résistants**, identifiés par un nom commençant par `R_` ou, après renommage, par une colonne commençant par `R` ;
- échantillons **sensibles**, identifiés par un nom commençant par `S_` ou par une colonne commençant par `S`.

Le pipeline suppose que les variants ont déjà été :

1. appelés individuellement, notamment avec `bcftools mpileup/call` d’après l’aide du script ;
2. annotés avec **SnpEff** ;
3. compressés en `VCF.gz` ;
4. indexés par un fichier `.tbi` ;
5. générés sur une même référence génomique.

Le pipeline réalise ensuite :

```text
VCF.gz individuels annotés SnpEff
        │
        ├── réétiquetage des noms d’échantillons
        ├── fusion multisample avec bcftools
        ├── conversion en TSV
        ├── calcul de l’AF à partir de AD et DP
        ├── filtres différentiels R/S
        ├── association aux gènes et exons du GTF
        ├── extraction des effets SnpEff
        └── visualisations :
                ├── fenêtres glissantes
                ├── variants individuels
                ├── agrégation par gène interactive
                └── agrégation par gène statique
```

## Portée exacte

Le dépôt **ne contient pas** les commandes complètes ayant produit les VCF individuels initiaux. L’appel de variants et l’annotation SnpEff sont des prérequis externes. Le code fourni commence réellement à partir de VCF individuels déjà annotés et indexés.

Deux méthodes différentes de rapprochement entre variants, gènes et effets sont présentes :

| Script | Fonction |
|---|---|
| `4_gtf_vcf_join.r` | jointure spatiale entre les positions des variants et les intervalles gène/exon d’un GTF simplifié |
| `4_variants_to_gene_effect_table.r` | extraction des annotations SnpEff `ANN`, association aux identifiants `HCON_*`, puis jointure avec les gènes du GTF |

Les deux scripts produisent des tables différentes et doivent être considérés comme deux sorties complémentaires ou deux approches concurrentes, et non comme une seule étape redondante parfaitement équivalente.

# Arborescence

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

Les dossiers numérotés `1_` à `7_` sont créés automatiquement au cours de l’exécution.

# Entrées et prérequis

## VCF individuels

Chaque sous-dossier de `0_samples/` doit contenir :

```text
<sample>.vcf.gz
<sample>.vcf.gz.tbi
```

Le script sélectionne le **premier fichier** correspondant à `*.vcf.gz` dans chaque dossier :

```bash
vcf="$(find "$sample" -maxdepth 1 -type f -name "*.vcf.gz" | head -n 1)"
```

Chaque dossier doit commencer par :

| Préfixe | Groupe |
|---|---|
| `R_` | résistant |
| `S_` | sensible |

Le nom du dossier devient le nouveau nom de l’échantillon dans le VCF.

## Référence GTF

Le premier fichier `*.gtf` trouvé directement dans `0_refseq/` est utilisé :

```bash
find ./0_refseq -maxdepth 1 -type f -name '*.gtf' | head -n 1
```

Le pipeline suppose notamment que :

- les chromosomes du GTF et du VCF utilisent les mêmes identifiants ;
- les attributs de gène suivent la syntaxe `gene_id "..."` ;
- les identifiants attendus peuvent suivre le motif `HCON_\d+`.

## Format des champs génotypiques

Le calcul de fréquence allélique suppose implicitement que les champs d’échantillon sont disposés comme suit après séparation par `:` :

```text
champ 1 : génotype ou autre champ
champ 2 : autre champ
champ 3 : DP
champ 4 : AD
```

Le code utilise :

```awk
dp = a[3]
ad = a[4]
```

L’ordre exact du champ `FORMAT` n’est toutefois pas lu dynamiquement. Le pipeline dépend donc d’un format de VCF spécifique.

# Logiciels et dépendances

Aucune version logicielle n’est explicitement fixée dans cette compilation.

| Outil ou package | Usage |
|---|---|
| Bash | orchestration générale |
| bcftools | lecture des noms d’échantillons, réétiquetage, indexation et fusion des VCF |
| SnpEff | annotation préalable des conséquences dans le champ `ANN` |
| awk | conversion, calcul d’AF, filtres et association initiale au GTF |
| gunzip | décompression du VCF fusionné |
| curl ou wget | téléchargement de fichiers `.tar` auxiliaires |
| R / Rscript | traitements tabulaires et figures |
| data.table | lecture rapide, jointures d’intervalles et agrégations |
| stringr | extraction des identifiants et annotations |
| readr | lecture et écriture TSV |
| dplyr | transformations et agrégations |
| ggplot2 | graphiques statiques |
| ggrepel | étiquetage des gènes |
| viridis | palette de densité |
| plotly | graphiques interactifs |
| htmlwidgets | export HTML |
| rmarkdown | vérification de Pandoc et export HTML autonome |
| Pandoc | génération des HTML `selfcontained` |

Le script télécharge également :

```text
bcftools.tar
r-plotly.tar
r-ggplot2.tar
```

depuis une URL IONOS. Cependant, ces archives ne sont **ni extraites ni ajoutées au `PATH` ou à la bibliothèque R** dans le script fourni. L’exécution repose donc en pratique sur des installations déjà accessibles dans l’environnement, sauf intervention externe non documentée.

# Orchestration générale

## Aide

```bash
bash VarEffPlot.sh -h
```

Options disponibles pour la figure statique par gène :

| Option | Effet |
|---|---|
| `-l <bp>` | borne génomique gauche |
| `-r <bp>` | borne génomique droite |
| `-b <valeur>` | borne inférieure de `AF_R - AF_S` |
| `-n` | ajout des identifiants de gènes |
| `-f` | force l’axe Y entre 0 et 1 |

## Exemple d’appel

```bash
bash VarEffPlot.sh -l 1000000 -r 5000000 -b 0.2 -n
```

## Logique de reprise

Chaque étape est protégée par un test sur l’existence du dossier de sortie :

```bash
if [ ! -d "4_join" ]; then
    ...
fi
```

Ainsi :

- si le dossier existe, l’étape est sautée ;
- pour recalculer une étape, le README recommande de supprimer le dossier correspondant puis de relancer ;
- l’existence du dossier, même vide ou incomplet, suffit à empêcher la réexécution.

# Étapes détaillées

## 0. Téléchargement de fichiers auxiliaires

```bash
mkdir -p "0_tools"
```

Fonction :

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

Fichiers demandés :

```bash
download_tool bcftools.tar
download_tool r-plotly.tar
download_tool r-ggplot2.tar
```

## 1. Réétiquetage des VCF individuels

**Dossier de sortie :** `1_vcf_reheadered/`

Pour chaque dossier d’échantillon :

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

Le script suppose implicitement que chaque VCF contient un seul nom d’échantillon. Si `bcftools query -l` renvoie plusieurs lignes, la table de renommage construite par `printf` ne sera pas correcte.

## 2. Fusion des VCF

**Dossier de sortie :** `2_vcf_merged/`

```bash
bcftools merge \
    -Oz \
    -o ./2_vcf_merged/all.vcf.gz \
    $(find ./1_vcf_reheadered/ -type f -name '*.vcf.gz' -print0 |
      xargs -0 printf '%s ')

bcftools index -t ./2_vcf_merged/all.vcf.gz
```

La fusion est effectuée sans option supplémentaire explicitant le comportement sur les allèles multiallé­liques, les sites absents ou les génotypes manquants. Les valeurs par défaut de la version installée de `bcftools` s’appliquent.

## 3A. Simplification du GTF

**Dossier de sortie :** `3_gtf_filtered/`  
**Fichier :** `gtf.tsv`

Seules les lignes `gene` et `exon` sont conservées.

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

Format produit :

```text
chromosome    gene_id    type    start    end
```

avec :

- `g` pour une ligne de gène ;
- `e` pour une ligne d’exon.

## 3B. Conversion et filtrage du VCF fusionné

**Dossier de sortie :** `3_vcf_filtered/`  
**Fichier final :** `vcf.tsv`

### 3B.1. Retrait des métadonnées `##`

```bash
gunzip -c "$cur" |
    awk 'BEGIN{OFS="\t"} /^##/ {next} {print}' \
    > 01.no_double_hash.tsv
```

La ligne `#CHROM` est conservée.

### 3B.2. Sélection des colonnes

Le script conserve :

```text
#CHROM
POS
REF
ALT
INFO
colonnes des échantillons
```

Les colonnes suivantes sont supprimées :

```text
ID
QUAL
FILTER
FORMAT
```

La suppression de `FORMAT` est importante : les champs d’échantillons sont ensuite interprétés selon un ordre fixe, sans conserver la description de cet ordre.

### 3B.3. Comptage des impacts SnpEff

Le nombre d’occurrences est calculé directement dans `INFO` :

```awk
high     = gsub(/\|HIGH\|/, "", tmp)
moderate = gsub(/\|MODERATE\|/, "", tmp)
modifier = gsub(/\|MODIFIER\|/, "", tmp)
low      = gsub(/\|LOW\|/, "", tmp)
```

Colonnes ajoutées :

```text
HIGH
MODERATE
MODIFIER
LOW
```

Ces valeurs représentent le nombre de sous-chaînes correspondantes dans `INFO`, pas nécessairement le nombre d’effets biologiquement indépendants.

### 3B.4. Calcul des fréquences alléliques

Pour chaque échantillon :

```awk
dp = a[3]
ad = a[4]
split(ad, adp, ",")
altad = (length(adp)>=2 ? adp[2] : adp[1])
AF = altad / dp
```

La valeur est formatée avec deux décimales :

```awk
sprintf("%.2f", AF)
```

Cas manquants ou invalides :

```text
./.  → 0.00
.    → 0.00
DP ≤ 0 → 0.00
AD absent ou invalide → 0.00
```

Une absence de donnée est donc convertie en fréquence allélique nulle, ce qui n’est pas équivalent biologiquement à une observation certaine de l’allèle de référence.

### 3B.5. Profondeur

Le minimum et le maximum de `DP` entre les échantillons sont calculés :

```text
minDP
maxDP
```

Les colonnes suivantes sont créées mais restent systématiquement à `NA` :

```text
minQUAL
maxQUAL
```

La colonne `QUAL` du VCF ayant été supprimée plus tôt, aucune qualité de variant n’est réellement résumée.

### 3B.6. Filtre des échantillons sensibles

Un variant est supprimé si sa fréquence allélique dépasse `0.1` dans **au moins un** échantillon dont le nom commence par `S` :

```awk
if($c!="NA" && ($c+0)>0.1){ bad=1 }
```

Le seuil est strict :

```text
AF > 0.1
```

Une valeur égale à `0.1` est conservée.

### 3B.7. Filtre des échantillons résistants

Paramètre :

```bash
MIN_R_COUNT=1
```

Un variant est conservé si :

```text
AF > 0.4
```

dans au moins `MIN_R_COUNT` colonnes commençant par `R`.

Dans l’état fourni :

```text
au moins un échantillon résistant avec AF > 0.4
```

Une valeur égale à `0.4` ne passe pas le filtre.

### 3B.8. Ajout préliminaire des gènes

Une étape `awk` ajoute une colonne `HCON_ID` en cherchant les gènes du GTF contenant la position du variant.

```text
CHROM    HCON_ID    POS    REF    ALT    INFO    ...
```

Plusieurs gènes sont concaténés avec une virgule. Les variants hors gène reçoivent `NA`.

Cette colonne est ajoutée au fichier final `3_vcf_filtered/vcf.tsv`, mais le script R de jointure spatiale réalise ensuite sa propre association variant–gène à partir du GTF simplifié.

## 4A. Jointure spatiale VCF–GTF

**Script :** `0_Rscripts/4_gtf_vcf_join.r`  
**Sortie :** `4_join/var.tsv`

Packages :

```r
library(data.table)
library(stringr)
```

### Lecture et contrôles

Colonnes VCF obligatoires :

```r
c("CHROM", "POS", "REF", "ALT", "INFO")
```

Le GTF simplifié reçoit les noms :

```r
c("chromosome", "gene_id", "type", "start", "end")
```

### Jointure aux gènes

Chaque variant est transformé en intervalle ponctuel :

```r
vcf_pos <- vcf[, .(
    row_id__,
    chromosome = CHROM,
    vstart = POS,
    vend = POS
)]
```

Jointure :

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

### Détection des exons

Une seconde jointure teste si la position appartient à un exon du même gène.

Code de région :

| Valeur | Interprétation |
|---|---|
| `e` | variant dans un exon |
| `g` | variant dans un gène, mais sans recouvrement exonique détecté |
| chaîne vide | variant non associé à un gène |

Le code ne distingue donc pas explicitement intron, UTR, promoteur ou intergénique dans `region_code`.

### Identifiants de gènes dans les annotations `MODIFIER`

La fonction `extract_modifier_gene_ids()` :

1. extrait `ANN=` depuis `INFO` ;
2. sépare les annotations sur les virgules ;
3. conserve uniquement celles dont le troisième champ vaut `MODIFIER` ;
4. cherche les motifs `HCON_\d+` dans les champs 4 et 5 ;
5. concatène les identifiants uniques avec `$`.

Sortie :

```text
modifier_gene_ids
```

### Colonnes de tête

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

Les autres colonnes du TSV VCF sont ensuite conservées.

## 4B. Table variant–gène–effet

**Script :** `0_Rscripts/4_variants_to_gene_effect_table.r`  
**Sortie principale :** `4_gene_var_eff/outnew`

Sorties intermédiaires :

```text
4_gene_var_eff/prepared.gtf.csv
4_gene_var_eff/prepared.vcf.csv
```

### Préparation du GTF

Seules les lignes `feature == "gene"` sont conservées.

Attributs extraits :

```text
gene_id
gene_name
```

Le préfixe `gene:` est supprimé du `gene_id`. Si `gene_name` est absent, le `gene_id` est utilisé.

### Extraction du champ SnpEff `ANN`

Le champ `ANN` est séparé selon les 16 champs standards attendus :

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

Le script extrait tous les motifs :

```r
HCON_\d+
```

depuis `gene_id` et `gene_name`.

Une ligne est créée par combinaison :

```text
variant × effet SnpEff × HCON_ID
```

Les identifiants qui ne figurent pas dans le GTF sont retirés.

### Table finale

Colonnes :

```text
chromosome
gene
gene_start
gene_end
variant_pos
variant_ref
variant_alt
<colonnes échantillons>
snpeff_effect
snpeff_impact
```

Cette sortie ne semble pas être utilisée par les scripts graphiques suivants, qui lisent principalement `4_join/var.tsv`.

## 5. Fenêtres glissantes et densité de variants

**Script :** `5_var_density.r`  
**Entrée :** `4_join/var.tsv`  
**Sortie :** `5_var_density/sliding_AF_R_minus_S_by_chr_plotly.html`

Packages :

```r
data.table
ggplot2
viridis
plotly
htmlwidgets
```

Paramètres effectifs :

```r
window_size <- 500000
step_size   <- 10000
```

Le commentaire en tête indique un pas de `100000 bp`, mais la variable exécutée utilise `10000 bp`.

### Détection des groupes

```r
s_cols <- grep("^S_", names(dt), value = TRUE)
r_cols <- grep("^R_", names(dt), value = TRUE)
```

Contrairement à certains autres scripts qui utilisent `^S` et `^R`, celui-ci exige spécifiquement des noms commençant par `S_` et `R_`.

### Calcul par variant

```r
mean_AF_S = moyenne des colonnes sensibles
mean_AF_R = moyenne des colonnes résistantes
delta_AF  = mean_AF_R - mean_AF_S
```

### Agrégation par fenêtre

Pour chaque chromosome :

```r
starts <- seq(1, max_pos, by = step_size)
win_end <- win_start + window_size - 1
```

Les fenêtres sont chevauchantes, car :

```text
fenêtre = 500 kb
pas = 10 kb
```

Mesures :

```text
n_variants
mean_delta
```

Les fenêtres vides reçoivent `n_variants = 0` et `mean_delta = NA`.

### Figure

- un panneau par chromosome ;
- courbe noire de `mean_delta` ;
- ligne horizontale à zéro ;
- trapèzes sous la courbe colorés par `n_variants` ;
- export Plotly autonome.

## 6. Visualisation interactive par variant

**Script :** `6_plotly_var.r`  
**Entrée :** `4_join/var.tsv`

Sorties :

```text
6_plotly_var/<chromosome>.html
6_plotly_var/top10pct_variants_by_absDiffAF.tsv
```

### Colonnes AF

```r
af_cols   <- grep("^[RS]", names(df), value = TRUE)
r_samples <- grep("^R", names(df), value = TRUE)
s_samples <- grep("^S", names(df), value = TRUE)
```

Cette règle peut sélectionner toute colonne commençant par `R` ou `S`, pas uniquement les échantillons, si d’autres noms de colonnes partagent ces initiales.

### Impact principal

Si les compteurs existent, la priorité est :

```text
HIGH > MODERATE > MODIFIER > LOW > NONE
```

### Différence de fréquence

```r
mean_R  <- moyenne des AF résistantes
mean_S  <- moyenne des AF sensibles
diff_AF <- mean_R - mean_S
```

### Top 10 %

Seuil :

```r
quantile(abs(diff_AF), probs = 0.90)
```

Les variants dont `|diff_AF|` est supérieur ou égal au 90e percentile sont exportés.

### Figure

Une bulle correspond à une ligne de la table et donc potentiellement à une combinaison variant–gène lorsque la jointure a produit plusieurs lignes pour un même variant.

Axes :

```text
X = position du variant
Y = moyenne(AF_R) - moyenne(AF_S)
```

Couleur :

```text
impact principal
```

Arrière-plan :

- rectangles bleus pour les intervalles de gènes ;
- rectangles verts pour les exons si leurs bornes existent ;
- sinon, traits ponctuels aux positions annotées `e`.

Un fichier HTML autonome est produit par chromosome. Pandoc est obligatoire.

## 7A. Visualisation statique agrégée par gène

**Script :** `7_ggplot_gene.r`  
**Entrée :** `4_join/var.tsv`

Sorties :

```text
7_ggplot_gene/ggplot_<chromosome>.pdf
7_ggplot_gene/top10pct_genes_by_absDiffAF.tsv
```

Seuls les variants avec :

```r
region_code %in% c("g", "e")
```

sont conservés.

### Agrégation

Groupement :

```text
chromosome
gene
gene_start
gene_end
```

Mesures :

```text
gene_mid
gene_len_bp
n_variants
AF moyenne par échantillon
best_impact
mean_R
mean_S
diff_AF
variants_per_kb
```

Densité :

```r
variants_per_kb = n_variants * 1000 / gene_len_bp
```

Le nombre de variants uniques repose sur :

```text
position:REF:ALT
```

### Top 10 % des gènes

```r
quantile(abs(gene_summary$diff_AF), probs = 0.90)
```

### Figure

```text
X = milieu du gène en Mb
Y = AF_R moyenne - AF_S moyenne
taille = variants par kb
couleur = meilleur impact
```

Dans l’état fourni, toutes les catégories d’impact utilisent toutefois la même couleur hexadécimale :

```r
"#440154"
```

Les commentaires à droite des couleurs suggèrent que d’autres couleurs avaient été envisagées, mais elles ne sont pas actives.

Les labels de gènes sont ajoutés avec `ggrepel` lorsque l’option `-n` est fournie.

## 7B. Visualisation interactive agrégée par gène

**Script :** `7_plotly_gene.r`  
**Entrée :** `4_join/var.tsv`

Sorties attendues dans :

```text
7_plotly_gene/
```

Ce script reprend l’agrégation par gène, mais ajoute une extraction plus détaillée de `INFO/ANN`.

### Sélection d’une annotation SnpEff

Pour chaque combinaison variant–gène, `extract_best_annotation_fast()` évalue les annotations selon :

1. correspondance du gène ;
2. correspondance de l’allèle alternatif ;
3. biotype `protein_coding` ;
4. présence d’un HGVS protéique ;
5. priorité de l’effet ;
6. niveau d’impact SnpEff.

Scores d’effet notamment :

| Effet | Score |
|---|---:|
| `stop_gained` | 100 |
| `frameshift_variant` | 95 |
| `start_lost` | 90 |
| site d’épissage canonique | 88 |
| insertion/délétion in-frame | 80 |
| `missense_variant` | 70 |
| synonyme | 60 |
| région d’épissage | 50 |
| UTR | 39–40 |
| intron | 30 |
| non-coding | 20 |
| upstream/downstream | 9–10 |
| intergénique | 5 |

Si un HGVS protéique est disponible, le préfixe `p.` est retiré et la notation devient le libellé d’effet affiché.

### Agrégation par gène

Le script calcule notamment :

```text
nombre de variants
longueur du gène
variants par kb
AF moyenne par échantillon
AF moyenne R
AF moyenne S
diff_AF
impact principal
liste des variants et effets retenus
```

Les figures interactives permettent l’inspection des informations par survol.

# Sorties principales

| Dossier | Fichier ou type de sortie |
|---|---|
| `1_vcf_reheadered/` | VCF individuels renommés et index `.tbi` |
| `2_vcf_merged/` | `all.vcf.gz`, `all.vcf.gz.tbi` |
| `3_gtf_filtered/` | `gtf.tsv` |
| `3_vcf_filtered/` | étapes intermédiaires `01` à `06`, puis `vcf.tsv` |
| `4_join/` | `var.tsv` |
| `4_gene_var_eff/` | `prepared.gtf.csv`, `prepared.vcf.csv`, `outnew` |
| `5_var_density/` | figure HTML de fenêtres glissantes |
| `6_plotly_var/` | HTML par chromosome et top 10 % des variants |
| `7_plotly_gene/` | figures HTML par chromosome |
| `7_ggplot_gene/` | PDF par chromosome et top 10 % des gènes |

# Seuils et paramètres critiques

| Paramètre | Valeur |
|---|---:|
| AF sensible maximale autorisée | aucun échantillon S avec `AF > 0.1` |
| AF résistante minimale | au moins un échantillon R avec `AF > 0.4` |
| `MIN_R_COUNT` | 1 |
| précision AF exportée | 2 décimales |
| fenêtre glissante | 500 000 bp |
| pas effectif | 10 000 bp |
| top variants | décile supérieur de `|AF_R - AF_S|` |
| top gènes | décile supérieur de `|AF_R - AF_S|` |

# Points de vigilance

1. **Versions absentes.** Aucune version de bcftools, SnpEff, R ou des packages R n’est consignée.
2. **Appel de variants externe.** Les commandes `mpileup/call`, leurs paramètres et les filtres antérieurs ne figurent pas dans le dépôt.
3. **Annotation SnpEff externe.** La base SnpEff et sa version ne sont pas indiquées.
4. **Archives téléchargées mais non installées.** Les `.tar` de bcftools et des packages R ne sont pas extraites ni utilisées explicitement.
5. **Premier fichier seulement.** Le premier GTF et le premier VCF de chaque dossier sont sélectionnés sans contrôle d’ambiguïté.
6. **Hypothèse d’un VCF mono-échantillon.** Le réétiquetage ne gère pas explicitement plusieurs noms dans un fichier source.
7. **FORMAT supprimé.** Le calcul de DP et AD repose sur les positions 3 et 4 des champs, pas sur les clés du champ `FORMAT`.
8. **Données manquantes transformées en AF = 0.** Cette convention peut renforcer artificiellement le contraste R/S.
9. **AF arrondie avant filtrage et analyses.** La précision est réduite à deux décimales.
10. **QUAL non exploitée.** `minQUAL` et `maxQUAL` restent `NA`.
11. **Aucun filtre explicite sur DP.** `minDP` et `maxDP` sont calculés mais aucun seuil de profondeur n’est appliqué.
12. **Filtres asymétriques.** Un seul R au-dessus de 0.4 suffit, alors qu’un seul S au-dessus de 0.1 exclut le variant.
13. **Préfixes variables.** Certains scripts cherchent `R_`/`S_`, d’autres simplement `R`/`S`.
14. **Jointure cartésienne possible.** Un variant chevauchant plusieurs gènes produit plusieurs lignes.
15. **Code `g` large.** Toute position dans un gène mais hors exon reçoit `g`, sans distinction intron/UTR.
16. **Deux méthodes d’association.** `4_join/var.tsv` et `4_gene_var_eff/outnew` ne sont pas strictement équivalents.
17. **Comptage des impacts par motif.** `HIGH`, `MODERATE`, etc. sont des occurrences de texte dans `INFO`.
18. **Fenêtres fortement chevauchantes.** 500 kb avec un pas de 10 kb induisent une forte autocorrélation visuelle.
19. **Incohérence commentaire/code.** Le commentaire annonce un pas de 100 kb, le code utilise 10 kb.
20. **Top 10 % fondé sur les lignes disponibles.** Les égalités au quantile peuvent produire plus de 10 % de lignes.
21. **Figures par variant potentiellement dupliquées.** Une jointure multi-gène peut afficher plusieurs bulles pour le même variant.
22. **Couleurs statiques identiques.** Le graphique ggplot par gène ne distingue pas visuellement les impacts dans l’état fourni.
23. **Reprise basée uniquement sur les dossiers.** Un dossier incomplet empêche la réexécution automatique.
24. **Pas de journalisation globale.** Les versions, dates, stdout/stderr et codes de sortie ne sont pas centralisés.
25. **Pas de contrôle des chromosomes compatibles.** Une différence de nomenclature entre GTF et VCF produit des variants non associés.

# Informations à ajouter pour le dépôt de publication

| Élément | Statut |
|---|---|
| version de bcftools | à renseigner |
| version de SnpEff et base utilisée | à renseigner |
| version de R et des packages | à renseigner si historiquement disponible |
| génome et annotation de référence | accession et version à renseigner |
| commandes d’appel de variants en amont | à documenter séparément |
| filtres appliqués avant ce pipeline | à documenter |
| définition biologique exacte des groupes R et S | à ajouter |
| tableau échantillon–condition | à ajouter |
| critères de profondeur et qualité | à préciser, ou indiquer explicitement leur absence |
| commit associé aux figures publiées | à fixer |
| date d’exécution | à ajouter |
| DOI et accession des données | à ajouter |

# Annexe — fichiers reproduits sans modification

Les fichiers ci-dessous sont reproduits tels qu’ils apparaissent dans la compilation. Les fichiers `.gitkeep` contenant uniquement l’indication `[Contenu binaire ignoré]` ne sont pas reproduits comme code analytique.


## `0_Rscripts/4_gtf_vcf_join.r`

Source brute : lignes 19–145 de la compilation.

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

Source brute : lignes 151–313 de la compilation.

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

Source brute : lignes 320–601 de la compilation.

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

Source brute : lignes 607–991 de la compilation.

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

Source brute : lignes 997–1313 de la compilation.

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

Source brute : lignes 1319–1742 de la compilation.

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

Source brute : lignes 1782–1794 de la compilation.

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

Source brute : lignes 1799–2140 de la compilation.

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
