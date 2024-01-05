


cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_WGS/mapping_eprinomectin_wgs


OUTPUT_DIR=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/ANALYSIS/WGS
find ~+ -type f -name '*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'


cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/ANALYSIS/WGS

## Eprinomectin resistance status
Bet - S
Chi - S
Luc - S
Ara - R
Bun - R
Mou - R


## within sample diversity
```bash

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

``

sed -i 's///g' hcontortus_eprinomectin_wgsfst.csv



Bet - S
Chi - S
Luc - S
Ara - R
Bun - R
Mou - R

chrom	start	end	snps	ARA.BET	ARA.BUN	ARA.CHI	ARA.LUC	ARA.MOU	BET.BUN	BET.CHI	BET.LUC	BET.MOU	BUN.CHI	BUN.LUC	BUN.MOU	CHI.LUC	CHI.MOU	LUC.MOU

ARA.BET
BET.BUN
BET.MOU
ARA.CHI
BUN.CHI
CHI.MOU
ARA.LUC
BUN.LUC
LUC.MOU



library(tidyverse)
rawdata <- read.table("hcontortus_eprinomectin_wgsfst.csv",  header=T)

data <- rawdata %>% pivot_longer(!chrom:snps, values_to="fst")

data_sus_v_res <- data %>% filter(str_detect(name, 'ARA.BET|BET.BUN|BET.MOU|ARA.CHI|BUN.CHI|CHI.MOU|ARA.LUC|BUN.LUC|LUC.MOU'))


ggplot(data_sus_v_res) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)")


data_sus_v_res_chr5 <- data_sus_v_res %>% filter(str_detect(chrom, 'hcontortus_chr5_Celeg_TT_arrow_pilon'))

ggplot(data_sus_v_res_chr5 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=37.5) + xlim(35,40)


data_sus_v_res_chr1 <- data_sus_v_res %>% filter(str_detect(chrom, 'hcontortus_chr1_Celeg_TT_arrow_pilon'))

ggplot(data_sus_v_res_chr1 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=7)


data_sus_v_res_chr4 <- data_sus_v_res %>% filter(str_detect(chrom, 'hcontortus_chr4_Celeg_TT_arrow_pilon'))

ggplot(data_sus_v_res_chr4 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)")



data_res_v_res <- data %>% filter(str_detect(name, 'ARA.BUN|ARA.MOU|BUN.MOU'))

data_res_v_res_chr5 <- data_res_v_res %>% filter(str_detect(chrom, 'hcontortus_chr5_Celeg_TT_arrow_pilon'))


ggplot(data_res_v_res) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)")


ggplot(data_res_v_res_chr5 ) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=37.5) + xlim(35,40)






data_sus_v_sus <- data %>% filter(str_detect(name, 'BET.CHI|BET.LUC|CHI.LUC'))

data_sus_v_sus_chr5 <- data_sus_v_sus %>% filter(str_detect(chrom, 'hcontortus_chr5_Celeg_TT_arrow_pilon'))


ggplot(data_sus_v_sus) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)")


ggplot(data_sus_v_sus_chr5) +
     geom_point(aes( start/1e6, fst, colour = chrom),  size = 0.1) +
     facet_grid(name ~ chrom , space="free_x", scales="free_x", switch="x") +
     theme_bw() +  theme(panel.spacing.x = unit(0, "lines"), legend.position = "none", text = element_text(size = 10), strip.text.y = element_text(size = 6)) +
     labs(x="Genomic position (Mb)", y = "Genetic differentiation (Fst)") +
     geom_vline(xintercept=37.5) + xlim(35,40)
