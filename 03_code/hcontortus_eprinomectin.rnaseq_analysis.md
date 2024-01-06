# hcontortus eprinomectin analysis: rna-seq

### author: Stephen Doyle

- RNA-seq analysis of eprinomectin-susceptible and resistant male and female parasites from six farms


## Data

```bash
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/RNASEQ

# use baton from Sam D to download fastq's from iRODs
module load ISG/singularity nextflow

bsub.py 10 baton "nextflow run /nfs/users/nfs_s/sd21/lustre_link/software/nextflow/baton-extractor/main.nf --study 7435 --runid 47866"

# run fastqc
bsub.py --queue long 10 fastqc "~sd21/bash_scripts/run_fastqc"

# once finished, run multiqc

mulltiqc .
```
[Multiqc report](../04_analysis/rnaseq_multiqc_report.html) 

| lane  | sample_ID |
|----------|---------|
| 48195#1  | CHI_M_1 |
| 48195#2  | CHI_M_2 |
| 48195#3  | CHI_M_3 |
| 48195#4  | CHI_F_1 |
| 48195#5  | CHI_F_2 |
| 48195#6  | CHI_F_3 |
| 48195#7  | BET_M_1 |
| 48195#8  | BET_M_2 |
| 48195#9  | BET_M_3 |
| 48195#10 | BET_F_1 |
| 48195#11 | BET_F_2 |
| 48195#12 | BET_F_3 |
| 48195#13 | LUC_M_1 |
| 48195#14 | LUC_M_2 |
| 48195#15 | LUC_M_3 |
| 48195#16 | LUC_F_1 |
| 48195#17 | LUC_F_2 |
| 48195#18 | LUC_F_3 |
| 48195#19 | ARA_M_1 |
| 48195#20 | ARA_M_2 |
| 48195#21 | ARA_M_3 |
| 48195#22 | ARA_F_1 |
| 48195#23 | ARA_F_2 |
| 48195#24 | ARA_F_3 |
| 48195#25 | BUN_M_1 |
| 48195#26 | BUN_M_2 |
| 48195#27 | BUN_M_3 |
| 48195#28 | BUN_F_1 |
| 48195#29 | BUN_F_2 |
| 48195#30 | BUN_F_3 |
| 48195#31 | MOU_M_1 |
| 48195#32 | MOU_M_2 |
| 48195#33 | MOU_M_3 |
| 48195#34 | MOU_F_1 |
| 48195#35 | MOU_F_2 |
| 48195#36 | MOU_F_3 |


## Reference
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/REF

# get genome
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.genomic.fa.gz

# get annotation
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.annotations.gff3.gz

# get transcripts
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS18.mRNA_transcripts.fa.gz

# unzip
for i in *gz; do 
    gunzip ${i}; 
    done

```




## Kalliso
```bash
cd /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ

# load kallisto (kallisto 0.50.1)
conda activate kallisto


# index the transcripts
ln -s /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/REF/haemonchus_contortus.PRJEB506.WBPS18.mRNA_transcripts.fa TRANSCRIPTS.fa

kallisto index --index TRANSCRIPTS.ixd TRANSCRIPTS.fa



DATA_DIR=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/RNASEQ


# run kallisto
DATA_DIR=/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/RAW/RNASEQ

while read LANE SAMPLE_ID; do \
    kallisto quant \
     --index TRANSCRIPTS.ixd \
     --output-dir kallisto_${SAMPLE_ID}_out \
     --bootstrap-samples 100 \
     --threads 7 \
     ${DATA_DIR}/${LANE}_1.fastq.gz ${DATA_DIR}/${LANE}_2.fastq.gz >>kallisto_${SAMPLE_ID}.log 2>&1;
done < ${DATA_DIR}/lane_sample_IDs.list &

```


## Sleuth
```bash

/nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/ANALYSIS/RNASEQ


```
- where "rnaseq_metadata.txt" contains:


| sample  | sex    | drug        | path                                                                                                                 |
|---------|--------|-------------|----------------------------------------------------------------------------------------------------------------------|
| ARA_F_1 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_ARA_F_1_out/abundance.h5 |
| ARA_F_2 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_ARA_F_2_out/abundance.h5 |
| ARA_F_3 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_ARA_F_3_out/abundance.h5 |
| ARA_M_1 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_ARA_M_1_out/abundance.h5 |
| ARA_M_2 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_ARA_M_2_out/abundance.h5 |
| ARA_M_3 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_ARA_M_3_out/abundance.h5 |
| BET_F_1 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BET_F_1_out/abundance.h5 |
| BET_F_2 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BET_F_2_out/abundance.h5 |
| BET_F_3 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BET_F_3_out/abundance.h5 |
| BET_M_1 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BET_M_1_out/abundance.h5 |
| BET_M_2 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BET_M_2_out/abundance.h5 |
| BET_M_3 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BET_M_3_out/abundance.h5 |
| BUN_F_1 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BUN_F_1_out/abundance.h5 |
| BUN_F_2 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BUN_F_2_out/abundance.h5 |
| BUN_F_3 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BUN_F_3_out/abundance.h5 |
| BUN_M_1 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BUN_M_1_out/abundance.h5 |
| BUN_M_2 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BUN_M_2_out/abundance.h5 |
| BUN_M_3 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_BUN_M_3_out/abundance.h5 |
| CHI_F_1 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_CHI_F_1_out/abundance.h5 |
| CHI_F_2 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_CHI_F_2_out/abundance.h5 |
| CHI_F_3 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_CHI_F_3_out/abundance.h5 |
| CHI_M_1 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_CHI_M_1_out/abundance.h5 |
| CHI_M_2 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_CHI_M_2_out/abundance.h5 |
| CHI_M_3 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_CHI_M_3_out/abundance.h5 |
| LUC_F_1 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_LUC_F_1_out/abundance.h5 |
| LUC_F_2 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_LUC_F_2_out/abundance.h5 |
| LUC_F_3 | female | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_LUC_F_3_out/abundance.h5 |
| LUC_M_1 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_LUC_M_1_out/abundance.h5 |
| LUC_M_2 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_LUC_M_2_out/abundance.h5 |
| LUC_M_3 | male   | susceptible | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_LUC_M_3_out/abundance.h5 |
| MOU_F_1 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_MOU_F_1_out/abundance.h5 |
| MOU_F_2 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_MOU_F_2_out/abundance.h5 |
| MOU_F_3 | female | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_MOU_F_3_out/abundance.h5 |
| MOU_M_1 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_MOU_M_1_out/abundance.h5 |
| MOU_M_2 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_MOU_M_2_out/abundance.h5 |
| MOU_M_3 | male   | resistant   | /nfs/users/nfs_s/sd21/lustre_link/haemonchus_contortus/EPRINOMECTIN/MAPPING_RNASEQ/kallisto_MOU_M_3_out/abundance.h5 |


