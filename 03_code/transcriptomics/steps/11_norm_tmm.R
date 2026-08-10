#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(optparse); library(edgeR)}))
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


cnts <- read.csv(op$`in`, row.names=1, header=TRUE)
GeneID <- rownames(cnts); Length <- cnts[,1]; cnt <- as.matrix(cnts[, -1])
TotalMapped <- colSums(cnt)
dge <- DGEList(counts = cnt)
dge <- calcNormFactors(dge, method = "TMM")
cnt <- cpm(dge, normalized.lib.sizes = TRUE)
colnames(cnt) <- colnames(cnt)
cnt <- cnt * TotalMapped / 1e6
out <- cbind(GeneID, Length, cnt)
write.csv(out, file="out/cnt.tmm.csv", row.names=FALSE, quote=FALSE)
