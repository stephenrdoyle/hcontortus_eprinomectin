#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(optparse); library(plotly)}))
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


dge <- read.csv(op$`in`, row.names=1)
if (!"log2FoldChange" %in% names(dge) || !"padj" %in% names(dge)) stop("DESeq2 columns not found")
dge$log2FC <- dge$log2FoldChange; dge$neglog10padj <- -log10(dge$padj)
dge <- subset(dge, neglog10padj > 0 & !is.na(neglog10padj))


volcan <- plot_ly() %>% layout(
title = "Volcano Plot",
xaxis = list(title = "Log2 Fold Change"),
yaxis = list(title = "-log10(Adjusted p-value)"),
showlegend = FALSE,
shapes = list(
list(type = "line", x0 = -1, y0 = 0, x1 = -1, y1 = max(dge$neglog10padj), line = list(dash = "dash")),
list(type = "line", x0 = 1, y0 = 0, x1 = 1, y1 = max(dge$neglog10padj), line = list(dash = "dash")),
list(type = "line", x0 = -5, y0 = 1.3, x1 = 5, y1 = 1.3, line = list(dash = "dash"))
)
)
threshold <- subset(dge, abs(log2FC) < 1 | neglog10padj < 1.3)
underexp <- subset(dge, log2FC <= -1 & neglog10padj >= 1.3)
overexp <- subset(dge, log2FC >= 1 & neglog10padj >= 1.3)
volcan <- volcan %>% add_trace(data = threshold, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
marker = list(size = 5), text = rownames(threshold))
volcan <- volcan %>% add_trace(data = underexp, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
marker = list(size = 5), text = rownames(underexp))
volcan <- volcan %>% add_trace(data = overexp, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
marker = list(size = 5), text = rownames(overexp))
htmlwidgets::saveWidget(volcan, "out/volcano_plot.html")
