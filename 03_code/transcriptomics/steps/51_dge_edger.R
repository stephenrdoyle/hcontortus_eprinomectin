#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(optparse); library(edgeR)}))
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


counts <- read.csv(op$`in`, row.names=1, header=TRUE)
counts <- counts[, -1] # remove Length
counts[] <- lapply(counts, function(x) as.numeric(x))


SeleCols <- function(data, pattern){ if (pattern!="") data[, grep(pattern, colnames(data)), drop=FALSE] else data }
ExCols <- function(data, pattern){ if (pattern!="") data[, -grep(pattern, colnames(data)), drop=FALSE] else data }


reference <- SeleCols(ExCols(counts, ""), "_S_")
comparison <- SeleCols(ExCols(counts, ""), "_R_")
X <- cbind(reference, comparison)


group <- factor(c(rep("Sensitive", ncol(reference)), rep("Resistant", ncol(comparison))), levels=c("Sensitive","Resistant"))
dge <- DGEList(counts=X, group=group)
design <- model.matrix(~group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit)
res <- topTags(qlf, n=Inf)
write.csv(as.data.frame(res), "out/dge.edger.csv")
