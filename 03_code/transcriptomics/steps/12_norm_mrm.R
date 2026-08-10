#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(optparse); library(DESeq2)}))
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


cnts <- read.csv(op$`in`, row.names=1, header=TRUE)
GeneID <- rownames(cnts); Length <- cnts[,1]; cnt <- as.matrix(cnts[, -1])
rownames(cnt) <- GeneID; colData <- data.frame(condition = colnames(cnt)); design <- ~ 1
dds <- DESeqDataSetFromMatrix(countData = cnt, colData = colData, design = design)
dds <- estimateSizeFactors(dds)
cnt <- counts(dds, normalized = TRUE)
out <- cbind(GeneID, Length, cnt)
write.csv(out, file="out/cnt.mrm.csv", row.names=FALSE, quote=FALSE)
