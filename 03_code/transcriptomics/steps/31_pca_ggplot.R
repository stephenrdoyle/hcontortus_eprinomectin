#!/usr/bin/env Rscript
suppressWarnings(suppressMessages({library(optparse); library(ggplot2)}))
opt_list <- list(make_option("--in", type="character"), make_option("--IdxCd", type="integer", default=1))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))


cnts <- read.csv(op$`in`, row.names=1, header=TRUE)
df <- t(as.matrix(cnts))
pca <- prcomp(df, scale=TRUE)
varexp <- (pca$sdev^2)/sum(pca$sdev^2)*100
split <- strsplit(rownames(pca$x), "_")
mat <- do.call(rbind, lapply(split, function(x){ if (length(x)>1) x else c(x, NA) }))
colnames(mat) <- paste0("cd", seq_len(ncol(mat)))
df_pca <- as.data.frame(cbind(pca$x, mat))


ncolname <- paste0("cd", op$IdxCd)
if (!(ncolname %in% names(df_pca))) ncolname <- names(df_pca)[grepl("^cd", names(df_pca))][1]


p <- ggplot(df_pca, aes(x = PC1, y = PC2, color = .data[[ncolname]])) +
  geom_point(size = 2) +
  coord_fixed() +
  labs(
    x = paste0("PC1 (", round(varexp[1], 2), "%)"),
    y = paste0("PC2 (", round(varexp[2], 2), "%)"),
    color = ncolname
  ) +
  theme_bw(base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "top"
  )


dir.create("out", recursive = TRUE, showWarnings = FALSE)
ggsave("out/pca.pdf", plot=p, width=6, height=6, device = cairo_pdf)
