# Récupération des options
args = commandArgs(trailingOnly=TRUE)
# Aide
if ("-h" %in% args | "-help" %in% args) {
  cat("🧬oooooooooooooo🧬", "\n")
  cat("🧬 MuGeCoExDeAn 🧬", "\n")
  cat("🧬°°°°°°°°°°°°°°🧬", "\n")
  cat("v.2025.03", "\n")
  cat("Robin LIOUTAUD, UMR1436-InTheRes-INRAE-Univ. Toulouse", "\n")
  cat("🆘▼▼▼ HELP MANUAL ▼▼▼", "\n")
  cat("🔣 Load and select gene count file(s)", "\n")
  cat("└ -i       : path of counts file, or folder in which folders in which count files, to load", "\n")
  cat("└ -k       : only keep samples with folder name containing that string", "\n")
  cat("└ -r       : remove all samples with folder name containing that string", "\n")
  cat("🧮 Counts normalization", "\n")
  cat("└ -tpm     : normalize in Transcripts (or Counts) Per Million", "\n")
  cat("└ -tmm     : normalize in Trimmed Mean of M-values (edgeR required)", "\n")
  cat("└ -mrm     : normalize in Mean of Ratio Method (DESeq2 required)", "\n")
  cat("📈 Post-normalization transformation", "\n")
  cat("└ -pntf    : méthode de transformation ap norm : rlog ou vst", "\n")
  cat("🔍 Exploration plot (PCA and clustered heatmap)", "\n")
  cat("└ -html ou -plotly : génère les PCA interactives multi-conditions", "\n")
  cat("└ -ggplot ou -pdf  : génère une PCA PC1 contre PC2 pour une condition donnée", "\n")
  cat("└ -IdxCd           : pour ggplot, rang de la condition à intituler et coloriser", "\n")
  cat("└ -chm ou -heatmap : génère la clustered heatmap", "\n")
  cat("📊 Plot counts by gene ID(s)", "\n")
  cat("└ -g      : LeGeneID", "\n")
  cat("└ -gl     : un fichier avec plusieurs GeneID en lignes", "\n")
  cat("└ -rank   : rang (ou indice) de la condition à représenter dans le sample name", "\n")
  cat("             exemple, pour l'échant. R.M.1, la condition M est en position 2", "\n")
  cat("└ -m      : afficher les médianes sur le graphique", "\n")
  cat("└ -s      : afficher les stats sur le graphique, * ** ***", "\n")
  cat("⚖️  Differential Gene Expression analysis", "\n")
  cat("└ -DESeq2 : DGE with DESeq2 (DESeq2 required. Will do MRM norm so don't do before.)", "\n")
  cat("└ -edgeR  : DGE with edgeR (edgeR required. TMM norm recommanded before.)", "\n")
  stop()
}

# Data loading ############################
t <- " Loading of gene count table(s) ..."
cat("⏳", t, "\r")
start <- Sys.time()
if ("-i" %in% args) {
  path <- args[which(args == "-i") + 1]
  cnts <- NULL
  Length <- NULL
  if (!file.info(path)$isdir) { # --- MODE FICHIER ---
    cnts <- read.csv(path, row.names = 1, header = TRUE)

    # ---------- CSV EXPORT : counts bruts combinés ----------
    try(write.csv(cnts, file = "./cnt.raw.csv", row.names = TRUE, quote = FALSE), silent = TRUE)

    # ---------- FILTRAGE COLONNES PAR -k / -r EN MODE FICHIER ----------
    cn <- colnames(cnts)
    keep <- rep(TRUE, length(cn))

    # Conserver toujours Length si présent (colonne 1 logique dans ton pipeline)
    if ("Length" %in% cn) {
      keep[which(cn == "Length")[1]] <- TRUE
    }

    # -k : ne garder que ce qui matche (peut être multiple, ET logique)
    for (i in seq_along(args)) {
      if (args[i] == "-k" && i + 1 <= length(args)) {
        pattern <- args[i + 1]
        keep <- keep & grepl(pattern, cn)
        # Ré-autoriser Length si écrasée par hasard
        if ("Length" %in% cn) keep[which(cn == "Length")[1]] <- TRUE
      }
    }
    # -r : exclure ce qui matche (peut être multiple, ET logique)
    for (i in seq_along(args)) {
      if (args[i] == "-r" && i + 1 <= length(args)) {
        pattern <- args[i + 1]
        keep <- keep & !grepl(pattern, cn)
        # Ré-autoriser Length si écrasée par hasard
        if ("Length" %in% cn) keep[which(cn == "Length")[1]] <- TRUE
      }
    }
    cnts <- cnts[, keep, drop = FALSE]

    # Re-export après filtrage (optionnel mais utile pour audit)
    try(write.csv(cnts, file = "./cnt.raw.filtered.csv", row.names = TRUE, quote = FALSE), silent = TRUE)

  } else { # --- MODE DOSSIER ---
    cntpaths <- list()
    samples <- sort(list.dirs(path, full.names = TRUE, recursive = FALSE))
    for (sample in samples) { # Récupération des path des counts
      cntfiles <- list.files(sample, pattern = "\\.cnt$", full.names = TRUE)
      if (length(cntfiles) > 0) {
        cntpaths <- c(cntpaths, cntfiles[which.min(nchar(basename(cntfiles)))])
      }
    }
    # Sélection/filtrage par -k
    for (arg in 1:length(args)) {
      if (args[arg] == "-k") {
        pattern <- args[arg + 1]
        cntpaths <- cntpaths[grepl(pattern, cntpaths)]
      }
    }
    # Exclusion par -r
    for (arg in 1:length(args)) {
      if (args[arg] == "-r") {
        pattern <- args[arg + 1]
        cntpaths <- cntpaths[!grepl(pattern, cntpaths)]
      }
    }
    cntslist <- list() # Récupération et listage des différents counts
    Length <- as.numeric(read.csv(cntpaths[[1]], header = TRUE)$Length)
    for (cntpath in cntpaths) {
      cnt <- read.csv(cntpath, header = TRUE)
      colnames(cnt)[1] <- "gene"
      cnt <- cnt[, !(names(cnt) %in% "Length")]
      colnames(cnt)[2] <- basename(dirname(cntpath))
      cntslist <- append(cntslist, list(cnt))
    }
    cnts <- cntslist[[1]]
    for (cnt in cntslist[-1]) {
      cnts <- merge(cnts, cnt, by = "gene", all = TRUE)
    }
    cnts <- cbind(cnts[, 1, drop = FALSE], Length, cnts[, -1, drop = FALSE]) # insère col Length en 2e position
    rownames(cnts) <- cnts[, 1]
    cnts <- cnts[, -1]
    cnts <- as.data.frame(cnts)

    # ---------- CSV EXPORT : counts bruts combinés (mode dossier) ----------
    try(write.csv(cnts, file = "./cnt.raw.csv", row.names = TRUE, quote = FALSE), silent = TRUE)
  }
}
stop <- Sys.time()
cat("✅", t, "completed in", stop - start, "seconds.", dim(cnts)[1], "x", dim(cnts)[2], "table.", "\n")

# Normalization #############################################
if ("-tpm" %in% args | "-tmm" %in% args | "-mrm" %in% args) {
  t <- " Normalizing gene counts table ..."
  cat("⏳", t, "\r")
  start <- Sys.time()
  GeneID <- rownames(cnts)
  Length <- cnts[, 1]
  cnt <- as.matrix(cnts[, -1])
  cnt[!sapply(cnt, is.numeric)] <- 0
  for (arg in args) {
    if (arg == "-tpm") {
      for (col in colnames(cnt)) {
        for (lin in 1:nrow(cnt)) {
          cnt[lin, col] <- cnt[lin, col] / Length[lin]
        }
        colsum <- sum(cnt[, col])
        for (lin in 1:nrow(cnt)) {
          cnt[lin, col] <- cnt[lin, col] / colsum * 1e6
        }
      }
    }
    if (arg == "-tmm") {
      suppressMessages(library(edgeR))
      TotalMapped <- colSums(cnt)
      dge <- DGEList(counts = cnt)
      dge <- calcNormFactors(dge, method = "TMM")
      cnt <- cpm(dge, normalized.lib.sizes = TRUE)
      colnames(cnt) <- colnames(cnt)
      cnt <- cnt * TotalMapped / 1e6
    }
    if (arg == "-mrm") { # Mean of Ratio Method (DESeq2)
      suppressMessages(library(DESeq2))
      rownames(cnt) <- GeneID
      colData <- data.frame(condition = colnames(cnt))
      design <- ~ 1
      dds <- DESeqDataSetFromMatrix(countData = cnt, colData = colData, design = design)
      dds <- estimateSizeFactors(dds)
      cnt <- counts(dds, normalized = TRUE)
    }
  }
  cnts <- cbind(GeneID, Length, cnt)
  cnts <- cnts[order(rownames(cnts)), ]
  write.csv(cnts, file = "./cnt.norm.csv", row.names = FALSE, quote = FALSE)
  stop <- Sys.time()
  cat("✅", t, "completed in", stop - start, "seconds.", "\n")
}

# Exploration ##################################################################################################################
if ("-html" %in% args | "-plotly" %in% args | "-ggplot" %in% args | "-pdf" %in% args | "chm" %in% args | "-heatmap" %in% args) {
  t <- " Exploring counts data ..."
  cat("⏳", t, "\r")
  filename <- "exploration"
  # df d'entrée pour PCA = cnt (matrix gènes x échantillons plus tard transposée)
  df <- cnt

  # ---------- CSV EXPORT : matrice entrée PCA (samples x genes) ----------
  try(write.csv(t(df), file = "./exploration.pca.matrix.csv", row.names = TRUE), silent = TRUE)

  # Fonction QC
  DfCheckBadData <- function(df) {
    cat("Total NA   :", sum(is.na(df)), "\n")
    cat("Total NULL :", sum(sapply(df, function(x) is.null(x))), "\n")
    cat("Total 0    :", sum(df == 0, na.rm = TRUE), "\n")
  }
  DfCheckBadData(df)

  df <- t(as.matrix(df))
  samples <- rownames(df)

  # Post-normalisation transformation
  if ("-pntf" %in% args) {
    pntf <- args[which(args == "-pntf")[1] + 1]
    if (pntf == "rlog") {
      epsilon <- 1
      df <- df + epsilon
      df <- log(df)
    } else if (pntf == "vst") {
      variances <- apply(df, 2, var)
      vst_transform <- sqrt(variances)
      df <- log(df + 1) / vst_transform
    }
    # ---------- CSV EXPORT : matrice transformée ----------
    try(write.csv(df, file = "./exploration.pca.transformed.csv", row.names = TRUE), silent = TRUE)
  }

  DfCheckBadData(df)
  # Retirer gènes constants/vides
  df <- df[, apply(df, 2, function(x) !all(is.na(x)) && length(unique(x)) > 1)]

  # PCA
  pca <- prcomp(df, scale = TRUE)
  pca$filename <- matrix(filename, nrow = 1, ncol = 1)

  # Conditions depuis les noms d'échantillons
  split_files <- strsplit(samples, "_")
  max_length <- max(sapply(split_files, length))
  conditions <- vector("list", max_length)
  for (i in seq_len(max_length)) {
    elements_at_i <- sapply(split_files, function(x) ifelse(length(x) >= i, x[i], NA))
    conditions[[i]] <- sort(unique(elements_at_i))
  }
  pca$conditions <- conditions
  pca$samples <- rownames(pca$x)
  pca$var <- pca$sdev^2
  pca$totvar <- sum(pca$var)
  pca$varexp <- list()
  pca$varexp$PC <- factor(1:length(pca$var))
  pca$varexp$pct <- (pca$var / pca$totvar) * 100
  pca$varexp$pctcum <- cumsum(pca$varexp$pct)

  # Table scores + conditions
  split <- strsplit(rownames(pca$x), "_")
  mat <- do.call(rbind, lapply(split, function(x) { if (length(x) > 1) return(x) else return(c(x, NA)) }))
  colnames(mat) <- paste0("cd", seq_len(ncol(mat)))
  pca$xc <- cbind(pca$x, mat)

  # ---------- CSV EXPORT : PCA ----------
  try(write.csv(pca$x, file = "./pca.scores.csv", row.names = TRUE), silent = TRUE)
  try(write.csv(pca$rotation, file = "./pca.loadings.csv", row.names = TRUE), silent = TRUE)
  try(write.csv(data.frame(PC = as.integer(pca$varexp$PC), pct = pca$varexp$pct, pctcum = pca$varexp$pctcum),
                file = "./pca.variance.csv", row.names = FALSE), silent = TRUE)
  try(write.csv(pca$xc, file = "./pca.scores_with_conditions.csv", row.names = TRUE), silent = TRUE)

  # Plotly (optionnel)
  if ("-html" %in% args | "-plotly" %in% args) {
    suppressMessages(library(plotly))
    suppressMessages(library(htmltools))
    plots <- list()

    plot_eigen <- plot_ly(
      data = as.data.frame(pca$var),
      x = pca$varexp$PC,
      y = pca$varexp$pct,
      type = "bar",
      name = "Explained Variance"
    ) %>% layout(
      xaxis = list(title = "Composantes principales"),
      yaxis = list(title = "% variance expliquée"),
      title = "Valeurs propres des PCs (scree plot)"
    )
    plots[[length(plots) + 1]] <- plot_eigen

    for (p in 1:((min(6, length(pca$sdev))) - 1)) {
      for (condition in pca$conditions) {
        cond <- paste(condition, collapse = "")
        plot_pca <- plot_ly() %>% layout(
          xaxis = list(title = list(text = paste0("PC", p, " (", round(pca$varexp$pct[p], 2), "%)"),
                                    font = list(size = 18, color = "black", weight = "bold"), standoff = 10)),
          yaxis = list(title = list(text = paste0("PC", p + 1, " (", round(pca$varexp$pct[p + 1], 2), "%)"),
                                    font = list(size = 18, color = "black", weight = "bold"))),
          legend = list(orientation = 'h', y = 1.18, x = 0.5, xanchor = "center",
                        font = list(size = 18, color = "black", weight = "bold"))
        )
        ExistingColors <- rainbow(length(condition))
        ExistingColors <- adjustcolor(ExistingColors, alpha.f = 1, red.f = 0.8, green.f = 0.8, blue.f = 0.8)
        i <- 1
        for (case in condition) {
          color <- ExistingColors[i]
          if (cond == 'FM') {
            if (case == 'F') color <- "coral"
            if (case == 'M') color <- "darkcyan"
          }
          if (cond == 'RS') {
            if (case == 'R') color <- "red"
            if (case == 'S') color <- "green"
          }
          pattern <- paste0("(^", case, "\\_)|(\\_", case, "$)|(\\_", case, "\\_)")
          pcasub <- pca$x[grepl(pattern, rownames(pca$x)), ]
          name <- case
          name <- gsub("\\bR\\b", "Resistant", name)
          name <- gsub("\\bS\\b", "Susceptible", name)
          name <- gsub("\\bM\\b", "Male", name)
          name <- gsub("\\bF\\b", "Female", name)
          plot_pca <- plot_pca %>% add_trace(
            x = pcasub[, paste0("PC", p)],
            y = pcasub[, paste0("PC", p + 1)],
            mode = "markers",
            marker = list(color = color, size = 7),
            text = rownames(pcasub),
            type = "scatter",
            textposition = "top",
            name = name
          )
          i <- i + 1
        }
        plot_pca <- plot_pca %>% config(toImageButtonOptions = list(format = 'svg', filename = 'pca'))
        plots[[length(plots) + 1]] <- plot_pca
      }
    }
    html_content <- tagList(div(style = "display: flex; flex-direction: column;",
                                lapply(plots, function(plot) div(plot))))
    save_html(html_content, file = paste0("./", pca$filename[1], ".pca.html"))
  }

  # ggplot PCA
  if ("-ggplot" %in% args | "-pdf" %in% args) {
    if ("-IdxCd" %in% args) {
      IdxCd = as.numeric(args[which(args == "-IdxCd")[1] + 1])
    }
    df_pca <- as.data.frame(pca$xc)
    suppressMessages(library(ggplot2))
    cond_cols <- grep("^cd", colnames(df_pca), ignore.case = TRUE)
    n <- IdxCd
    if (n <= length(cond_cols)) {
      ncolname <- cond_cols[IdxCd]
      uniqueconds <- unique(df_pca[[ncolname]])
      uniquecolors <- rainbow(length(uniqueconds))
      uniquecolors <- adjustcolor(uniquecolors, alpha.f = 1, red.f = 0.8, green.f = 0.8, blue.f = 0.8)
      for (uniquecond in 1:length(uniqueconds)) {
        df_pca$color[df_pca[[ncolname]] == uniqueconds[uniquecond]] <- uniquecolors[uniquecond]
      }
      df_pca$color[df_pca[[ncolname]] == "R"] <- "red"
      df_pca[[ncolname]][df_pca[[ncolname]] == "R"] <- "Resistants"
      df_pca$color[df_pca[[ncolname]] == "S"] <- "darkgreen"
      df_pca[[ncolname]][df_pca[[ncolname]] == "S"] <- "Susceptibles"
      df_pca$color[df_pca[[ncolname]] == "M"] <- "darkcyan"
      df_pca[[ncolname]][df_pca[[ncolname]] == "M"] <- "Males"
      df_pca$color[df_pca[[ncolname]] == "F"] <- "coral"
      df_pca[[ncolname]][df_pca[[ncolname]] == "F"] <- "Females"
    }
    df_pca$PC1 <- as.numeric(df_pca$PC1)
    df_pca$PC2 <- as.numeric(df_pca$PC2)
    xside <- max(abs(df_pca$PC1)); xside <- xside + xside/5
    yside <- max(abs(df_pca$PC2)); yside <- yside + yside/5
    p <- ggplot(df_pca, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = df_pca$color), size = 2) +
      coord_fixed(ratio=1) +
      scale_x_continuous(limits = c(-xside, xside), expand=c(0, 0)) +
      scale_y_continuous(limits = c(-yside, yside), expand=c(0, 0)) +
      labs(x = paste0("PC1 (", round(pca$varexp$pct[1], 2), "%)"),
           y = paste0("PC2 (", round(pca$varexp$pct[2], 2), "%)")) +
      scale_color_identity(
        breaks = unique(df_pca$color),
        labels = unique(df_pca[[paste0("cd", n)]]),
        guide = "legend"
      ) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "gray90", color = NA),
        plot.background = element_rect(fill = NA, color = NA),
        text = element_text(family = "sans", face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        aspect.ratio = 1
      )
    ggsave(filename = "./pca.pdf", plot = p, width = 6, height = 6, device = cairo_pdf)
  }

  # Heatmap / correlogram
  if ("-chm" %in% args | "-heatmap" %in% args | "-correlogram" %in% args) {
    suppressMessages(library(plotly))
    suppressMessages(library(pheatmap))
    suppressMessages(library(htmltools))

    # Renommage échantillons (clean labels)
    rownames(df)[rownames(df) == "R.ARA.M.1"] <- "ARA-M1"
    rownames(df)[rownames(df) == "R.ARA.M.2"] <- "ARA-M2"
    rownames(df)[rownames(df) == "R.ARA.M.3"] <- "ARA-M3"
    rownames(df)[rownames(df) == "R.ARA.F.1"] <- "ARA-F1"
    rownames(df)[rownames(df) == "R.ARA.F.2"] <- "ARA-F2"
    rownames(df)[rownames(df) == "R.ARA.F.3"] <- "ARA-F3"
    rownames(df)[rownames(df) == "R.BET.M.1"] <- "BET-M1"
    rownames(df)[rownames(df) == "R.BET.M.2"] <- "BET-M2"
    rownames(df)[rownames(df) == "R.BET.M.3"] <- "BET-M3"
    rownames(df)[rownames(df) == "R.BET.F.1"] <- "BET-F1"
    rownames(df)[rownames(df) == "R.BET.F.2"] <- "BET-F2"
    rownames(df)[rownames(df) == "R.BET.F.3"] <- "BET-F3"
    rownames(df)[rownames(df) == "R.BUN.M.1"] <- "BUN-M1"
    rownames(df)[rownames(df) == "R.BUN.M.2"] <- "BUN-M2"
    rownames(df)[rownames(df) == "R.BUN.M.3"] <- "BUN-M3"
    rownames(df)[rownames(df) == "R.BUN.F.1"] <- "BUN-F1"
    rownames(df)[rownames(df) == "R.BUN.F.2"] <- "BUN-F2"
    rownames(df)[rownames(df) == "R.BUN.F.3"] <- "BUN-F3"
    rownames(df)[rownames(df) == "R.MOU.M.1"] <- "MOU-M1"
    rownames(df)[rownames(df) == "R.MOU.M.2"] <- "MOU-M2"
    rownames(df)[rownames(df) == "R.MOU.M.3"] <- "MOU-M3"
    rownames(df)[rownames(df) == "R.MOU.F.1"] <- "MOU-F1"
    rownames(df)[rownames(df) == "R.MOU.F.2"] <- "MOU-F2"
    rownames(df)[rownames(df) == "R.MOU.F.3"] <- "MOU-F3"
    rownames(df)[rownames(df) == "S.CHI.M.1"] <- "CHI-M1"
    rownames(df)[rownames(df) == "S.CHI.M.2"] <- "CHI-M2"
    rownames(df)[rownames(df) == "S.CHI.M.3"] <- "CHI-M3"
    rownames(df)[rownames(df) == "S.CHI.F.1"] <- "CHI-F1"
    rownames(df)[rownames(df) == "S.CHI.F.2"] <- "CHI-F2"
    rownames(df)[rownames(df) == "S.CHI.F.3"] <- "CHI-F3"
    rownames(df)[rownames(df) == "S.LUC.M.1"] <- "LUC-M1"
    rownames(df)[rownames(df) == "S.LUC.M.2"] <- "LUC-M2"
    rownames(df)[rownames(df) == "S.LUC.M.3"] <- "LUC-M3"
    rownames(df)[rownames(df) == "S.LUC.F.1"] <- "LUC-F1"
    rownames(df)[rownames(df) == "S.LUC.F.2"] <- "LUC-F2"
    rownames(df)[rownames(df) == "S.LUC.F.3"] <- "LUC-F3"

    # Matrices de QC
    sampleCor <- cor(t(df))
    sampleDistsBrut <- 1 - sampleCor
    sampleDists <- as.dist(sampleDistsBrut)
    sampleDistMatrix <- as.matrix(sampleDists)

    # ---------- CSV EXPORT : QC ----------
    try(write.csv(sampleCor, file = "./exploration.sampleCor.csv", row.names = TRUE), silent = TRUE)
    try(write.csv(sampleDistsBrut, file = "./exploration.sampleDist.csv", row.names = TRUE), silent = TRUE)

    valmin <- min(sampleDistMatrix)
    valmax <- max(max(sampleDistMatrix), 0)

    palette <- c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")
    props <- seq(0, 1, length.out = length(palette))
    allcolors <- c()
    for (i in 1:(length(palette) - 1)) {
      num_colors_segment <- round((props[i + 1] - props[i]) * 100)
      palette_segment <- colorRampPalette(c(palette[i], palette[i + 1]))(num_colors_segment)
      allcolors <- c(allcolors, palette_segment)
    }
    breaks <- valmin + props * (valmax - valmin)
    legend_breaks <- breaks

    n <- max(nrow(sampleDistMatrix), ncol(sampleDistMatrix))
    fontsize <- floor(80 / sqrt(n))
    cell_size <- fontsize + floor(fontsize / 8)
    pdf_side <- floor((4 * cell_size + n * cell_size + fontsize * 8) / 60)

    pdf(paste0(filename, ".clustheatmap.pdf"), width = pdf_side, height = pdf_side, family = "sans")
    result <- pheatmap(
      sampleDistMatrix,
      clustering_distance_rows = sampleDists,
      clustering_distance_cols = sampleDists,
      border_color = NA,
      color = allcolors,
      breaks = seq(from = valmin, to = valmax, length.out = length(allcolors) + 1),
      legend_breaks = legend_breaks,
      legend_labels = round(legend_breaks, 3),
      fontsize = fontsize,
      fontsize_row = fontsize,
      fontsize_col = fontsize,
      fontsize_number = floor(fontsize * 2 / 3),
      cellwidth = cell_size,
      cellheight = cell_size
    )

    # ---------- CSV EXPORT : matrice de distances ordonnée par clustering ----------
    try({
      row_order <- labels(result$tree_row)
      col_order <- labels(result$tree_col)
      if (!is.null(row_order) && !is.null(col_order)) {
        clustered <- sampleDistMatrix[row_order, col_order, drop = FALSE]
        write.csv(clustered, file = "./exploration.sampleDist.clustered.csv", row.names = TRUE)
      }
    }, silent = TRUE)

    # Analyse outlier niveau 1
    row_clusters <- cutree(result$tree_row, k = 2)
    population_cluster <- names(row_clusters[row_clusters == 1])
    outlier_cluster <- names(row_clusters[row_clusters == 2])
    population <- as.data.frame(t(df))
    outlier <- population[, outlier_cluster, drop = FALSE]
    population <- population[, population_cluster, drop = FALSE]
    mean_population <- rowMeans(population, na.rm = TRUE)
    mean_outlier <- rowMeans(outlier, na.rm = TRUE)
    outlier_data <- numeric(length(mean_outlier))
    names(outlier_data) <- rownames(population)
    for (i in seq_along(mean_outlier)) {
      if (!is.na(mean_outlier[i]) && !is.na(mean_population[i]) && mean_outlier[i] > 0 && mean_population[i] > 0) {
        if (mean_outlier[i] >= mean_population[i]) {
          outlier_data[i] <- (mean_outlier[i] + 1) / (mean_population[i] + 1)
        } else {
          outlier_data[i] <- -(mean_population[i] + 1) / (mean_outlier[i] + 1)
        }
      } else {
        outlier_data[i] <- 0
      }
    }
    if (length(outlier_data) != nrow(population)) {
      stop("Problème dans le calcul : le nombre de gènes en sortie diffère de celui en entrée.")
    }
    top_genes <- sort(outlier_data, decreasing = TRUE)
    write.csv(data.frame(Gene = names(top_genes), Difference = top_genes),
              file = "./ClustNiv1GeneExprDiff.csv", row.names = FALSE)
    dev.off()
  }
  cat("✅", t, " completed.", "\n")
}

# Plot counts by gene ID ##########################################################
if ("-g" %in% args | "-gene" %in% args | "-gl" %in% args | "-genelist" %in% args) {
  t <- " Plotting counts for specific gene(s) ..."
  cat("⏳", t, "\r")
  start <- Sys.time()
  suppressMessages(library(dplyr))
  suppressMessages(library(ggplot2))

  colnames <- colnames(cnts)
  colnames <- colnames[-1]
  colnames <- as.character(colnames)
  colnames <- strsplit(colnames, split = "\\_")

  # FONCTION PAR GENE
  countplot <- function(cnt, args) {
    cnt <- cnt[rownames(cnt) == gene, , drop = FALSE]
    cnt <- t(cnt)
    counts <- cnt[, 1]
    theconditions <- as.data.frame(do.call(rbind, strsplit(rownames(cnt), "_")))
    for (j in seq_along(args)) { if (args[j] == "-rank" | args[j] == "-rang") { rang <- as.numeric(args[j+1]) } }
    if (is.null(rang)) { print("ERREUR : -r manquant") }
    condition <- as.factor(theconditions[, rang])
    samples <- theconditions[, apply(theconditions, 2, function(x) length(unique(x)) > 1)]
    samples <- apply(samples, 1, function(x) paste(x, collapse = "-"))
    row <- data.frame(Count = counts, Condition = condition)
    rownames(row) <- samples

    # ---------- CSV EXPORT : données du plot par gène ----------
    try(write.csv(row, file = paste0("./", gene, ".counts.by.condition.csv"), row.names = TRUE), silent = TRUE)

    # Couleurs / stats
    unicdt <- unique(row$Condition)
    paires <- combn(unicdt, 2, simplify = FALSE)
    colors <- setNames(rainbow(length(unicdt), s = 0.9, v = 0.7), unicdt)
    row$Color <- colors[as.character(row$Condition)]
    row$Color[row$Condition == "R"] <- "red"
    row$Color[row$Condition == "S"] <- "darkgreen"
    row$Color[row$Condition == "M"] <- "darkcyan"
    row$Color[row$Condition == "F"] <- "coral"

    p <- ggplot(row, aes(x = Condition, y = Count, color = Color)) + scale_x_discrete()
    if ("-m" %in% args | "-median" %in% args) {
      medians <- row %>% group_by(Condition) %>% summarise(median_value = median(Count))
      p <- p + geom_segment(data = medians,
                            aes(x = as.numeric(Condition) - 0.2, xend = as.numeric(Condition) + 0.2,
                                y = median_value, yend = median_value),
                            color = "black", linewidth = 1)
    }
    xtitle <- gene
    ytitle <- "Raw counts"
    norms <- list()
    for (arg in args) {
      if (arg == "-tpm" | arg == "-tmm" | arg == "-mrm" | arg == "rlog") {
        norms <- append(norms, gsub("-", "", arg))
      }
    }
    if (length(norms) > 0) {
      norms <- paste(norms, collapse = "->")
      ytitle <- paste0("Normalized counts (", norms, ")")
    }
    p <- p + geom_point(position = position_jitterdodge(jitter.width = 1), size = 3) +
      scale_color_identity() +
      theme_minimal() +
      labs(title = xtitle, x = "Condition", y = ytitle) +
      theme(
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.text = element_text(size = 15, color = "black", face = "bold"),
        axis.title = element_text(size = 15, color = "black", face = "bold"),
        plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5)
      )
    if ("-s" %in% args | "-stat" %in% args) {
      library(ggsignif)
      test_results <- lapply(paires, function(pair) {
        data_subset <- row[row$Condition %in% pair, ]
        test <- wilcox.test(Count ~ Condition, data = data_subset, exact = FALSE)
        return(list(pair = pair, p_value = test$p.value))
      })
      significant_pairs <- lapply(test_results, function(res) {
        if (res$p_value < 0.05) return(res$pair) else return(NULL)
      })
      significant_pairs <- significant_pairs[!sapply(significant_pairs, is.null)]
      p <- p + geom_signif(
        comparisons = significant_pairs,
        map_signif_level = TRUE,
        y_position = seq(from = max(row$Count)+1, by = 1*max(row$Count)/10, length.out = length(significant_pairs)),
        tip_length = 0.02,
        color = "black"
      )
    }
    suppressMessages(ggsave(paste0(gene, ".pdf"), plot = p))
  }

  # Exécution
  for (j in seq_along(args)) {
    if (args[j] == "-g" | args[j] == "-gene") {
      gene <- args[j+1]
      countplot(cnt, args)
    }
  }
  for (j in seq_along(args)) {
    if (args[j] == "-gl" | args[j] == "-genelist") {
      genelist <- readLines(args[j+1])
      genelist <- genelist[genelist != ""]
      if (length(genelist) > 0 && genelist[length(genelist)] == "") {
        genelist <- genelist[-length(genelist)]
      }
      genelist <- as.character(genelist)
      for (thegene in genelist) {
        gene <- thegene
        countplot(cnt, args)
      }
    }
  }
  stop <- Sys.time()
  cat("✅", t, "completed in", stop - start, "seconds.", "\n")
}

# Differential Gene Expression ################
if ("-DESeq2" %in% args | "-edgeR" %in% args) {
  t <- " Analyzing Differential Gene Expression ..."
  cat("⏳", t, "\r")
  start <- Sys.time()
  # DESeq2
  if ("-DESeq2" %in% args) {
    if ("-mrm" %in% args) {
      print("❌ don't use -mrm option before DESeq2 analysis,")
      print("   that will internally run MRM normalization.")
      stop()
    }
    suppressMessages(library(DESeq2))
    suppressMessages(library(ggplot2))
    filename <- "dgedeseq2"

    counts <- cnts
    counts <- counts[, -1] # supprime Length
    counts <- as.matrix(counts)
    counts <- matrix(as.integer(counts), nrow = nrow(counts), ncol = ncol(counts),
                     dimnames = list(rownames(counts), colnames(counts)))

    SeleCols <- function(data, pattern) {
      if (pattern != "") {
        cols <- grep(pattern, colnames(data))
      } else {
        cols <- 1:ncol(data)
      }
      return(data[, cols])
    }
    ExCols <- function(data, pattern) {
      if (pattern != "") {
        exclude_cols <- grep(pattern, colnames(data))
        data <- data[, -exclude_cols]
      }
      return(data)
    }

    reference <- counts
    reference <- ExCols(reference, "")
    reference <- SeleCols(reference, "_CHI_|_LUC_")

    comparison <- counts
    comparison <- ExCols(comparison, "")
    comparison <- SeleCols(comparison, "_ARA_|_BUN_|_MOU_")

    CountsRefComp <- cbind(reference, comparison)
    conditions <- factor(c(rep("ref", ncol(reference)), rep("comp", ncol(comparison))), levels=c("ref", "comp"))
    dds <- DESeqDataSetFromMatrix(countData = CountsRefComp, colData = data.frame(condition = conditions), design = ~ condition)

    # Normalisation & export
    dds <- estimateSizeFactors(dds)
    normcnt <- counts(dds, normalized=TRUE)
    write.csv(normcnt, file = paste0(filename, ".mrm.csv"), row.names = TRUE)

    # DGE
    dds <- DESeq(dds)
    res <- results(dds)
    res <- res[order(res$padj), ]
    write.csv(res, file = paste0(filename, ".deseq2.dge.csv"), row.names = TRUE)

    # Plots qualité
    pdf(paste0(filename, "_", "dispersion", ".pdf")); plotDispEsts(dds); dev.off()
    svg(paste0(filename, "_", "dispersion", ".svg"))
    dispersions <- dispersionFunction(dds)
    plot(dispersions, type = "b", col = "blue", pch = 16, cex = 1.5, main = "Custom Dispersion Plot")
    dev.off()

    # Nettoyage mrm intermédiaires à côté du path (si utilisés)
    mrm <- list.files(path, full.names = TRUE)
    mrm <- mrm[grep(".mrm.csv", mrm)]
    file.remove(mrm)
  }

  # edgeR
  if ("-edgeR" %in% args) {
    suppressMessages(library(edgeR))
    counts <- cnts
    counts <- counts[, -1]
    counts <- counts[, -1]
    counts <- as.data.frame(counts)
    counts[] <- lapply(counts, function(x) as.numeric(x))

    SeleCols <- function(data, pattern) {
      if (pattern != "") {
        cols <- grep(pattern, colnames(data))
      } else {
        cols <- 1:ncol(data)
      }
      return(data[, cols])
    }
    ExCols <- function(data, pattern) {
      if (pattern != "") {
        exclude_cols <- grep(pattern, colnames(data))
        data <- data[, -exclude_cols]
      }
      return(data)
    }

    counts <- ExCols(counts, "")
    counts <- SeleCols(counts, "")
    reference <- counts; reference <- ExCols(reference, ""); reference <- SeleCols(reference, "_S_")
    comparison <- counts; comparison <- ExCols(comparison, ""); comparison <- SeleCols(comparison, "_R_")

    counts <- cbind(reference, comparison)
    group <- factor(c(rep("Sensitive", ncol(reference)), rep("Resistant", ncol(comparison))), levels = c("Sensitive", "Resistant"))
    dge <- DGEList(counts=counts, group=group)
    design <- model.matrix(~group)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    qlf <- glmQLFTest(fit)
    results <- topTags(qlf, n=Inf)
    write.csv(as.data.frame(results), "./resultats_DGE.csv")
  }
  stop <- Sys.time()
  cat("✅", t, "completed.", "\n")
}

# Volcano (optionnel) ##########################
if ("-volcano" %in% args) {
  library(ggplot2); library(tools); library(plotly)
  filename <- "volcano"
  dge <- res; dge <- as.data.frame(dge)
  names(dge)[names(dge) == "log2FoldChange"] <- "log2FC"
  names(dge)[names(dge) == "padj"] <- "neglog10padj"
  dge$neglog10padj <- -log10(dge$neglog10padj)
  dge <- subset(dge, neglog10padj > 0)
  dge <- dge[!is.na(dge$neglog10padj), ]
  volcan <- plot_ly() %>% layout(
    title = "Volcano Plot",
    xaxis = list(title = "Log2 Fold Change",
                 titlefont = list(color = "black", family = "serif", size = 22, weight = "bold"),
                 tickfont = list(color = "black", family = "serif", size = 22, weight = "bold")),
    yaxis = list(title = "-log10(Adjusted p-value)",
                 titlefont = list(color = "black", family = "serif", size = 22, weight = "bold"),
                 tickfont = list(color = "black", family = "serif", size = 22, weight = "bold")),
    showlegend = FALSE,
    shapes = list(
      list(type = "line", x0 = -1, y0 = 0, x1 = -1, y1 = max(dge$neglog10padj), line = list(color = "green", dash = "dash")),
      list(type = "line", x0 =  1, y0 = 0, x1 =  1, y1 = max(dge$neglog10padj), line = list(color = "green", dash = "dash")),
      list(type = "line", x0 = -5, y0 = 1.3, x1 = 5,  y1 = 1.3, line = list(color = "green", dash = "dash"))
    )
  )
  threshold <- subset(dge, abs(log2FC) < 1 | abs(neglog10padj) < 1.3)
  underexp <- subset(dge, log2FC <= -1 & neglog10padj >= 1.3)
  overexp  <- subset(dge, log2FC >=  1 & neglog10padj >= 1.3)
  volcan <- volcan %>% add_trace(data = threshold, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
                                 marker = list(size = 5, color = "black"), text = rownames(threshold))
  volcan <- volcan %>% add_trace(data = underexp, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
                                 marker = list(size = 5, color = "blue"), text = rownames(underexp))
  volcan <- volcan %>% add_trace(data = overexp, x = ~log2FC, y = ~neglog10padj, type = "scatter", mode = "markers",
                                 marker = list(size = 5, color = "red"), text = rownames(overexp))
  volcan <- volcan %>% config(toImageButtonOptions = list(format = 'svg', filename = 'volcano'))
  htmlwidgets::saveWidget(volcan, paste0(filename, "_", "volcano_plot.html"))
}

