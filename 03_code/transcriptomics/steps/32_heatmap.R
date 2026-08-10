#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(optparse); library(pheatmap)}))
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


cnts <- read.csv(op$`in`, row.names=1, header=TRUE)
df <- t(as.matrix(cnts[, -1]))
# sample distances
sampleCor <- cor(df)
sampleDistsBrut <- 1 - sampleCor
sampleDists <- as.dist(sampleDistsBrut)
sampleDistMatrix <- as.matrix(sampleDists)
write.csv(sampleCor, file = "out/exploration.sampleCor.csv")
write.csv(sampleDistsBrut, file = "out/exploration.sampleDist.csv")


valmin <- min(sampleDistMatrix); valmax <- max(max(sampleDistMatrix), 0)
pdf_side <- 8
pdf("out/exploration.clustheatmap.pdf", width=pdf_side, height=pdf_side)
res <- pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
border_color=NA)
# export clustered order if available
try({
row_order <- labels(res$tree_row); col_order <- labels(res$tree_col)
if (!is.null(row_order) && !is.null(col_order)) {
clustered <- sampleDistMatrix[row_order, col_order, drop=FALSE]
write.csv(clustered, file = "out/exploration.sampleDist.clustered.csv")
}
}, silent=TRUE)
dev.off()
