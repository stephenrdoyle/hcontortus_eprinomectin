


cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_WGS/mapping_eprinomectin_wgs


OUTPUT_DIR=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/ANALYSIS/WGS
find ~+ -type f -name '*.bam*' -exec ln -vs "{}" $OUTPUT_DIR/ ';'


cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/ANALYSIS/WGS




module load grenedalf/0.2.0

## within sample diversity
```bash
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



``



## between sample diversity 
```bash
grenedalf fst

grenedalf fst --write-sample-coverage --write-sample-alt-freq --write-total-frequency --sam-min-map-qual 30 --sam-min-base-qual 30 --file-prefix MHCO3_v_MHCO18_pools --separator-char tab --sam-path MHCO3_P0_L3_n200_01.bam --sam-path MHCO18_P0_L3_n200_IVM_01.bam


```
