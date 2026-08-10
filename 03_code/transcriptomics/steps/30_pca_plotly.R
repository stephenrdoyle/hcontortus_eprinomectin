#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
  library(plotly)
  library(htmltools)
}))

## --- PARSING DES OPTIONS ----
opt_list <- list(make_option("--in", type="character"))
op <- parse_args(OptionParser(option_list=opt_list))
stopifnot(file.exists(op$`in`))

dir.create("out", showWarnings = FALSE)

## --- LOAD ---
cnts <- read.csv(
  op$`in`,
  row.names   = 1,
  header      = TRUE,
  check.names = FALSE
)

# on enlève la première colonne (souvent "Length")
df <- t(as.matrix(cnts[, -1, drop = FALSE]))  # samples x genes

# garder une trace de la matrice utilisée
write.csv(
  t(df),
  file      = "out/exploration.pca.matrix.csv",
  row.names = TRUE
)

DfCheckBadData <- function(df){
  cat("Total NA:", sum(is.na(df)),
      " Total 0:", sum(df == 0, na.rm = TRUE), "\n")
}
DfCheckBadData(t(df))

df <- as.matrix(df)
samples <- rownames(df)

## --- FILTRAGE gènes à variance nulle ---
if (nrow(df) < 2) {
  stop("PCA: besoin d'au moins 2 samples (lignes) pour calculer une PCA.")
}

sd_per_gene <- apply(df, 2, function(x) sd(x, na.rm = TRUE))
keep <- !is.na(sd_per_gene) & sd_per_gene > 0

if (!any(keep)) {
  stop("PCA: toutes les colonnes (gènes) ont une variance nulle. Impossible de faire une PCA.")
}

if (sum(keep) < ncol(df)) {
  message("PCA: filtrage de ", sum(!keep),
          " gènes constants / variance nulle avant PCA (",
          sum(keep), " gènes conservés).")
  df <- df[, keep, drop = FALSE]
}

## --- PCA ---
pca <- prcomp(df, scale. = TRUE)
pca$var    <- pca$sdev^2
pca$totvar <- sum(pca$var)
varexp     <- data.frame(
  PC  = 1:length(pca$var),
  pct = (pca$var / pca$totvar) * 100
)

# sauvegardes standard
write.csv(pca$x,         file = "out/pca.scores.csv")
write.csv(pca$rotation,  file = "out/pca.loadings.csv")
write.csv(
  transform(varexp, pctcum = cumsum(pct)),
  file      = "out/pca.variance.csv",
  row.names = FALSE
)

## --- EXTRACTION DES CONDITIONS À PARTIR DES NOMS DE SAMPLES ---

# rownames(pca$x) = noms de samples
split <- strsplit(rownames(pca$x), "_")

mat <- do.call(
  rbind,
  lapply(split, function(x){
    if (length(x) > 1) x else c(x, NA)
  })
)
colnames(mat) <- paste0("cd", seq_len(ncol(mat)))

# data.frame scores + conditions
scores <- as.data.frame(pca$x)
scores$sample <- rownames(scores)
conds  <- as.data.frame(mat, stringsAsFactors = FALSE)
scores <- cbind(scores, conds)

# sauvegarde scores + conditions
write.csv(
  scores,
  file      = "out/pca.scores_with_conditions.csv",
  row.names = FALSE
)

## --- PLOTS PLOTLY ---

plots <- list()

# 1) Scree plot (inchangé)
plot_eigen <- plot_ly(
  data = varexp,
  x    = ~PC,
  y    = ~pct,
  type = "bar",
  name = "Explained Variance"
) %>%
  layout(
    xaxis = list(title = "PCs"),
    yaxis = list(title = "% variance"),
    title = "Scree plot"
  )

plots[[length(plots) + 1]] <- plot_eigen

# 2) PCA colored by each condition (cd1, cd2, ...)
maxpc    <- min(6, ncol(pca$x))  # max PC index affiché
pc_names <- colnames(pca$x)

# liste des colonnes conditions (cd1, cd2, ...)
cond_cols <- colnames(conds)

for (cond in cond_cols) {

  cond_vec <- as.factor(conds[[cond]])

  # si toute la condition est NA ou constante, on saute
  if (all(is.na(cond_vec)) || length(unique(na.omit(cond_vec))) <= 1) {
    next
  }

  for (p in 1:(maxpc - 1)) {

    pcx <- pc_names[p]
    pcy <- pc_names[p + 1]

    xvals <- scores[[pcx]]
    yvals <- scores[[pcy]]

    plot_pca <- plot_ly(
      x     = xvals,
      y     = yvals,
      type  = "scatter",
      mode  = "markers",
      text  = scores$sample,
      color = cond_vec,
      marker = list(size = 7)
    ) %>%
      layout(
        # pas de titre global pour la figure
        xaxis  = list(
          title     = paste0(
            pcx, " (", round(varexp$pct[p], 2), "%)"
          ),
          titlefont = list(size = 18),
          tickfont  = list(size = 14)
        ),
        yaxis  = list(
          title     = paste0(
            pcy, " (", round(varexp$pct[p + 1], 2), "%)"
          ),
          titlefont = list(size = 18),
          tickfont  = list(size = 14)
        ),
        legend = list(
          orientation = "h",
          y           = 1.12,
          x           = 0.5,
          xanchor     = "center",
          font        = list(size = 16)
        )
      ) %>%
      config(toImageButtonOptions = list(format = "svg", filename = "pca"))

    plots[[length(plots) + 1]] <- plot_pca
  }
}

# 3) Export HTML
html_content <- tagList(
  div(
    style = "display:flex;flex-direction:column;",
    lapply(plots, function(pl) div(pl))
  )
)

save_html(html_content, file = "out/exploration.pca.html")

