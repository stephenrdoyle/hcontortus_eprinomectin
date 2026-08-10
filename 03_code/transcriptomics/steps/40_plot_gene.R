#!/usr/bin/env Rscript
# passthrough of original flags to keep logic; expects -g or -gl and -rank
args <- commandArgs(trailingOnly=TRUE)
stopifnot(any(args %in% c("-g","-gene","-gl","-genelist")))


# decide counts input (prefer most recent normalized if exists)
choose_counts <- function(){
cands <- c("out/cnt.mrm.csv","out/cnt.tmm.csv","out/cnt.tpm.csv","out/cnt.raw.filtered.csv","out/cnt.raw.csv")
for (f in cands) if (file.exists(f)) return(f)
stop("No counts file found in ./out")
}
file <- choose_counts()
cnts <- read.csv(file, row.names=1, header=TRUE)
Length <- cnts[,1]; cnt <- as.matrix(cnts[, -1])


suppressMessages(library(ggplot2))


countplot <- function(gene, rank, with_median=FALSE, with_stats=FALSE){
sub <- t(cnt[rownames(cnt)==gene, , drop=FALSE])
counts <- sub[,1]
parts <- as.data.frame(do.call(rbind, strsplit(rownames(sub), "_")))
cond <- as.factor(parts[, rank])
samples <- apply(parts, 1, function(x) paste(x, collapse = "-"))
df <- data.frame(Count=counts, Condition=cond, Sample=samples)
write.csv(df, file = file.path("out", paste0(gene, ".counts.by.condition.csv")), row.names=TRUE)
p <- ggplot(df, aes(x=Condition, y=Count)) + geom_point(position=position_jitter(width=0.2), size=2) +
theme_minimal() + labs(title=gene, x="Condition", y="Counts")
if (with_median) {
med <- aggregate(Count~Condition, df, median)
p <- p + geom_segment(data=med, aes(x=as.numeric(Condition)-0.2, xend=as.numeric(Condition)+0.2,
y=Count, yend=Count), inherit.aes=FALSE)
}
ggsave(file.path("out", paste0(gene, ".pdf")), plot=p, width=6, height=4)
}


rank <- {
idx <- which(args %in% c("-rank","-rang")); if (length(idx)) as.integer(args[idx+1]) else stop("-rank missing")
}


if (any(args %in% c("-g","-gene"))) {
idx <- which(args %in% c("-g","-gene"))[1]; gene <- args[idx+1]
countplot(gene, rank, with_median = any(args %in% c("-m","-median")), with_stats = any(args %in% c("-s","-stat")))
}
if (any(args %in% c("-gl","-genelist"))) {
idx <- which(args %in% c("-gl","-genelist"))[1]; gl <- readLines(args[idx+1])
gl <- gl[gl!=""]
for (gene in gl) countplot(gene, rank, with_median = any(args %in% c("-m","-median")), with_stats = any(args %in% c("-s","-stat")))
}
