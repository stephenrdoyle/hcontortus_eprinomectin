#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library(optparse)))
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


cnts <- read.csv(op$`in`, row.names=1, header=TRUE)
GeneID <- rownames(cnts); Length <- cnts[,1]; cnt <- as.matrix(cnts[, -1])
cnt[!sapply(cnt, is.numeric)] <- 0
for (col in colnames(cnt)) {
cnt[, col] <- cnt[, col] / Length
colsum <- sum(cnt[, col])
cnt[, col] <- cnt[, col] / colsum * 1e6
}
out <- cbind(GeneID, Length, cnt)
write.csv(out, file="out/cnt.tpm.csv", row.names=FALSE, quote=FALSE)
