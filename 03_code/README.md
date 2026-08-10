# Structuration of the codes here

## Samples

Samples are designed in .manifest files.

## Genomics

The first bloc is about genomics.

genomics.sh is the analysis pipeline per-sample

genomics_VariantComparison is the script to compare allele frequency of variants in genes between the susceptible and the resistant groups

## Transcriptomics

The second bloc is about transcriptomics.

transcriptomics.sh is the analysis pipeline per-sample

MuGeCoExDeAn is the R tool created to generae DGE data and explore them

Upset plot on DGE data

## Warning

Some codes designed to be launched in specific environments : HPC cluster Sanger with LSF, R locally, ...
Need to be adapted for reproducibility in your environment.

## Authors

Robin Lioutaud, Stephen R. Doyle
