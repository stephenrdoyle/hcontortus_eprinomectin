#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library(optparse)))
opt_list <- list(make_option("--in", type="character"), make_option("--method", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


cnts <- read.csv(op$`in`, row.names=1, header=TRUE)
Length <- cnts[,1]
df <- t(as.matrix(cnts[, -1]))
if (op$method=="rlog") { df <- log(df + 1) } else if (op$method=="vst") {
variances <- apply(df, 2, var); df <- log(df + 1) / sqrt(variances)
} else stop("--method must be rlog or vst")
# export transformed sample x gene matrix for PCA
write.csv(df, file = "out/exploration.pca.transformed.csv", row.names = TRUE)
base <- sub(".*/|\\.csv$","", op$`in`)
outfile <- file.path("out", sprintf("%s.%s.csv", sub("\\.csv$","", base), op$method))
out <- cbind(GeneID=rownames(cnts), Length, t(df))
write.csv(out, file=outfile, row.names=FALSE, quote=FALSE)
