# hcontortus eprinomectin analysis: whole genome analyses

### author: Stephen Doyle

- pool-seq whole genome analysis of pools of parasites from two eprinomectin-susceptible and four resistant farms
    - Chi - Susceptible
    - Luc - Susceptible
    - Ara - Resistant
    - Bet - Resistant
    - Bun - Resistant
    - Mou - Resistant
- note that originally, there was meant to be three resistant and three susceptible farms
     - after sequencing and some analysis, we found that Bet looked more resistant. This was confirmed by some new phenotyping performed between sending the samples for sequencing and analysis. Hence, 4 resistant and 2 susceptible.



## Reference genome
- using the H. contortus reference genome hosted on WormBase Parasite 
     - (https://parasite.wormbase.org/Haemonchus_contortus_prjeb506/Info/Index/)
     - genome paper: https://doi.org/10.1038/s42003-020-01377-3 

- I used a local version already on my HPC, but the same genome can be downloaded as follows

```bash
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/REF

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa.gz
```


## Mapping
- mapping was performed using our in-house nextflow pipeline 
- requires a sample manifest, here called "eprinomectin_wgs.mapping.manifest", which is used in the mapping pipline. 
- It contains the following information, and links the sample name to the sequnecing data

| ID     | R1                                                                                           | R2                                                                                           |
|--------|----------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
| CHI_L3 | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#1_1.fastq.gz | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#1_2.fastq.gz |
| BET_L3 | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#2_1.fastq.gz | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#2_2.fastq.gz |
| LUC_L3 | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#3_1.fastq.gz | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#3_2.fastq.gz |
| ARA_L3 | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#4_1.fastq.gz | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#4_2.fastq.gz |
| BUN_L3 | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#5_1.fastq.gz | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#5_2.fastq.gz |
| MOU_L3 | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#6_1.fastq.gz | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/48115_2#6_2.fastq.gz |


```bash
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_WGS

# run mapping pipeline
bsub.py 10 mapping_wgs "mapping-helminth --reference HAEM_V4_final.chr.fa --input eprinomectin_wgs.mapping.manifest --outdir mapping_eprinomectin_wgs"


# collect and check mapping stats using multiqc
multiqc .

```
[Multiqc report](../04_analysis/wgs_multiqc_report.html) 






## Collect mapped data for analysis
```bash
# go to mapping directory
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_WGS/mapping_eprinomectin_wgs

# nice way to find all the bams, and them make symbolic links to a target directory
OUTPUT_DIR=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/ANALYSIS/WGS
find ~+ -type f -name '*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'
```



## within sample diversity
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/ANALYSIS/WGS


module load grenedalf/0.2.0

bsub.py --queue long 20 grenedlf_diversity \
"grenedalf diversity \
--file-prefix hcontortus_eprinomectin_wgs \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 50 \
--filter-sample-max-coverage 5000 \
--pool-sizes 1000 \
--window-type sliding \
--window-sliding-width 10000 \
--separator-char tab \
--sam-path ARA_L3.bam \
--sam-path BET_L3.bam \
--sam-path BUN_L3.bam \
--sam-path CHI_L3.bam \
--sam-path LUC_L3.bam \
--sam-path MOU_L3.bam"

sed -i 's/_L3:1//g' hcontortus_eprinomectin_wgsdiversity.csv
```

```R
library(tidyverse)
rawdata <- read.table("hcontortus_eprinomectin_wgsdiversity.csv",  header=T)

data <- rawdata %>% pivot_longer(!chrom:end, values_to="data")

ARA.snp_count	ARA.coverage_fraction	ARA.theta_pi_abs	ARA.theta_pi_rel	ARA.theta_watterson_abs	ARA.theta_watterson_rel	ARA.tajimas_d

data_snp_count <- data %>% filter(str_detect(name, 'snp_count'))
data_coverage_fraction <- data %>% filter(str_detect(name, 'coverage_fraction'))
data_theta_pi_abs <- data %>% filter(str_detect(name, 'theta_pi_abs'))
data_theta_pi_rel <- data %>% filter(str_detect(name, 'theta_pi_rel'))
data_theta_watterson_abs <- data %>% filter(str_detect(name, 'theta_watterson_abs'))
data_theta_watterson_rel <- data %>% filter(str_detect(name, 'theta_watterson_rel'))
data_tajimas_d <- data %>% filter(str_detect(name, 'tajimas_d'))



# data_snp_count
ggplot(data_snp_count, aes(start+5000/1e6, data,colour = chrom)) + 
    geom_point(size = 0.1) + 
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6))

# data_coverage_fraction
ggplot(data_coverage_fraction, aes(start+5000/1e6, data,colour = chrom)) + 
    geom_point(size = 0.1) + 
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6))

# data_theta_pi_abs
ggplot(data_theta_pi_abs, aes(start+5000/1e6, log10(data),colour = chrom)) + 
    geom_point(size = 0.1) + 
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6))

# data_theta_pi_rel
ggplot(data_theta_pi_rel, aes(start+5000/1e6, log10(data),colour = chrom)) + 
    geom_point(size = 0.1) + 
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6))

# data_theta_watterson_abs
ggplot(data_theta_watterson_abs, aes(start+5000/1e6, data,colour = chrom)) + 
    geom_point(size = 0.1) + 
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6))

ggplot(data_tajimas_d, aes(start+5000/1e6, data,colour = chrom)) + 
    geom_point(size = 0.1) + 
    facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
    theme_bw() + 
    theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6))


```





## between sample diversity
```bash

module load grenedalf/0.2.0

bsub.py --queue long 20 grenedlf_fst \
"grenedalf fst \
--method unbiased-nei \
--file-prefix hcontortus_eprinomectin_wgs \
--sam-min-map-qual 30 \
--sam-min-base-qual 30 \
--filter-sample-min-count 2 \
--filter-sample-min-coverage 50 \
--filter-sample-max-coverage 5000 \
--pool-sizes 1000 \
--window-type sliding \
--window-sliding-width 10000 \
--separator-char tab \
--sam-path ARA_L3.bam \
--sam-path BET_L3.bam \
--sam-path BUN_L3.bam \
--sam-path CHI_L3.bam \
--sam-path LUC_L3.bam \
--sam-path MOU_L3.bam"

sed -i 's/_L3:1//g' hcontortus_eprinomectin_wgsfst.csv
``

    - Chi - S
    - Luc - S
    - Ara - R
    - Bet - R
    - Bun - R
    - Mou - R


ARA.BET	ARA.BUN	ARA.CHI	ARA.LUC	ARA.MOU	BET.BUN	BET.CHI	BET.LUC	BET.MOU	BUN.CHI	BUN.LUC	BUN.MOU	CHI.LUC	CHI.MOU	LUC.MOU

```R
library(tidyverse)
rawdata <- read.table("hcontortus_eprinomectin_wgsfst.csv",  header=T)

data <- rawdata %>% pivot_longer(!chrom:snps, values_to="fst")

data_sus_v_res <- data %>% filter(str_detect(name, 'ARA.CHI|BET.CHI|BUN.CHI|CHI.MOU|ARA.LUC|BUN.LUC|LUC.MOU|BET.LUC'))


ggplot(data_sus_v_res) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     ylim(0,1)


# zoom on chromosome 1
data_sus_v_res_chr1 <- data_sus_v_res %>% filter(str_detect(chrom, 'hcontortus_chr1_Celeg_TT_arrow_pilon'))

ggplot(data_sus_v_res_chr1 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=7)

# zoom on chromosome 2
data_sus_v_res_chr2 <- data_sus_v_res %>% filter(str_detect(chrom, 'hcontortus_chr2_Celeg_TT_arrow_pilon'))

ggplot(data_sus_v_res_chr2 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)")


# zoom on chromosome 5
data_sus_v_res_chr4 <- data_sus_v_res %>% filter(str_detect(chrom, 'hcontortus_chr4_Celeg_TT_arrow_pilon'))

ggplot(data_sus_v_res_chr4 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)")

# zoom on chromosome 5
data_sus_v_res_chr5 <- data_sus_v_res %>% filter(str_detect(chrom, 'hcontortus_chr5_Celeg_TT_arrow_pilon'))

ggplot(data_sus_v_res_chr5 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)")

ggplot(data_sus_v_res_chr5 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=37.5) + xlim(35,40)






# resistant vs resistant

data_res_v_res <- data %>% filter(str_detect(name, 'ARA.BET|ARA.BUN|ARA.MOU|BUN.MOU|BET.BUN|BET.MOU'))

data_res_v_res_chr5 <- data_res_v_res %>% filter(str_detect(chrom, 'hcontortus_chr5_Celeg_TT_arrow_pilon'))


ggplot(data_res_v_res) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     ylim(0,1)


ggplot(data_res_v_res_chr5 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=37.5) + xlim(35,40)




# susceptible vs susceptible

data_sus_v_sus <- data %>% filter(str_detect(name, 'CHI.LUC'))

data_sus_v_sus_chr5 <- data_sus_v_sus %>% filter(str_detect(chrom, 'hcontortus_chr5_Celeg_TT_arrow_pilon'))


ggplot(data_sus_v_sus) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     ylim(0,1)


ggplot(data_sus_v_sus_chr5) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=37.5) + xlim(35,40)


```