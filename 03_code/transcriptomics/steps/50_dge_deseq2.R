#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(optparse); library(DESeq2); library(ggplot2)}))
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


counts <- read.csv(op$`in`, row.names=1, header=TRUE)
counts <- counts[, -1] # remove Length
counts <- as.matrix(counts)
counts <- matrix(as.integer(counts), nrow=nrow(counts), ncol=ncol(counts), dimnames=list(rownames(counts), colnames(counts)))


SeleCols <- function(data, pattern){ if (pattern!="") data[, grep(pattern, colnames(data)), drop=FALSE] else data }
ExCols <- function(data, pattern){ if (pattern!="") data[, -grep(pattern, colnames(data)), drop=FALSE] else data }


reference <- SeleCols(ExCols(counts, "_F\\."), "_CTL_")
comparison <- SeleCols(ExCols(counts, "_F\\."), "_MOX_")
CountsRefComp <- cbind(reference, comparison)
conditions <- factor(c(rep("ref", ncol(reference)), rep("comp", ncol(comparison))), levels=c("ref","comp"))
dds <- DESeqDataSetFromMatrix(countData = CountsRefComp, colData = data.frame(condition=conditions), design = ~ condition)


# normalization export
dds <- estimateSizeFactors(dds)
normcnt <- counts(dds, normalized=TRUE)
write.csv(normcnt, file="out/dgedeseq2.mrm.csv")


# DGE
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
write.csv(res, file="out/dge.deseq2.dge.csv")


pdf("out/dgedeseq2_dispersion.pdf"); plotDispEsts(dds); dev.off()
