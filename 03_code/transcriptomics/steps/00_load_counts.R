#!/usr/bin/env Rscript
suppressWarnings(suppressMessages(library(optparse)))
opt_list <- list(
make_option("--path", type="character"),
make_option("--keep", type="character", default=NA),
make_option("--remove", type="character", default=NA)
)
op <- parse_args(OptionParser(option_list=opt_list))


dir.create("out", showWarnings = FALSE)


msg <- function(...) cat("[load]", ..., "\n")
stopifnot(!is.null(op$path))


if (!file.info(op$path)$isdir) {
cnts <- read.csv(op$path, row.names=1, header=TRUE)
try(write.csv(cnts, file="out/cnt.raw.csv", row.names=TRUE, quote=FALSE), silent=TRUE)
cn <- colnames(cnts); keep <- rep(TRUE, length(cn))
if ("Length" %in% cn) keep[which(cn=="Length")[1]] <- TRUE
if (!is.na(op$keep)) {
keep <- keep & grepl(op$keep, cn)
if ("Length" %in% cn) keep[which(cn=="Length")[1]] <- TRUE
}
if (!is.na(op$remove)) {
keep <- keep & !grepl(op$remove, cn)
if ("Length" %in% cn) keep[which(cn=="Length")[1]] <- TRUE
}
cnts <- cnts[, keep, drop=FALSE]
try(write.csv(cnts, file="out/cnt.raw.filtered.csv", row.names=TRUE, quote=FALSE), silent=TRUE)
} else {
cntpaths <- list(); samples <- sort(list.dirs(op$path, full.names=TRUE, recursive=FALSE))
for (sample in samples) {
cntfiles <- list.files(sample, pattern = "\\.cnt$", full.names = TRUE)
if (length(cntfiles) > 0) cntpaths <- c(cntpaths, cntfiles[which.min(nchar(basename(cntfiles)))])
}
if (!is.na(op$keep)) cntpaths <- cntpaths[grepl(op$keep, cntpaths)]
if (!is.na(op$remove)) cntpaths <- cntpaths[!grepl(op$remove, cntpaths)]
cntslist <- list(); Length <- as.numeric(read.csv(cntpaths[[1]], header=TRUE)$Length)
for (cntpath in cntpaths) {
cnt <- read.csv(cntpath, header=TRUE)
colnames(cnt)[1] <- "gene"; cnt <- cnt[, !(names(cnt) %in% "Length")]
colnames(cnt)[2] <- basename(dirname(cntpath))
cntslist <- append(cntslist, list(cnt))
}
cnts <- cntslist[[1]]
for (cnt in cntslist[-1]) cnts <- merge(cnts, cnt, by="gene", all=TRUE)
cnts <- cbind(cnts[,1,drop=FALSE], Length, cnts[, -1, drop=FALSE])
rownames(cnts) <- cnts[,1]; cnts <- cnts[, -1]; cnts <- as.data.frame(cnts)
try(write.csv(cnts, file="out/cnt.raw.csv", row.names=TRUE, quote=FALSE), silent=TRUE)
}
msg("done")
