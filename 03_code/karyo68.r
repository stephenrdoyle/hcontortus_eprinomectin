#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  library(grid)
  library(scales)
  library(dplyr)
  # ggh4x est optionnel (axes Y différents par rangée)
  # install.packages("ggh4x") si besoin
})

# =========================
# CLI
# =========================
option_list <- list(
  make_option(c("--dge"),     type = "character", help = "DGE CSV brut"),
  make_option(c("--fai"),     type = "character", help = "FAI (genomic.fa.fai)"),
  make_option(c("--gtf"),     type = "character", help = "GTF (geneset.gtf)"),
  make_option(c("--regions"), type = "character", default = NA, help = "CSV régions: Chr,startMb,endMb (sans header)"),
  make_option(c("--title"),   type = "character", default = " ", help = "Titre plot"),
  make_option(c("--out"),     type = "character", default = "karyoplot.pdf", help = "PDF de sortie"),
  make_option(c("--padj"),    type = "double",  default = 0.05, help = "Seuil FDR"),
  make_option(c("--lfc"),     type = "double",  default = 0.58, help = "Seuil |log2FC|"),
  make_option(c("--win"),     type = "integer", default = 1e6,  help = "Fenêtre hypergéom (bp)"),
  make_option(c("--step"),    type = "integer", default = 1e5,  help = "Pas (bp)"),
  make_option(c("--width"),   type = "double",  default = 16,   help = "Largeur PDF (in)"),
  make_option(c("--height"),  type = "double",  default = 5,    help = "Hauteur PDF (in)"),
  make_option(c("--hyperHighlight"), type = "logical", default = TRUE,
              help = "TRUE/FALSE: colorer aussi les régions significatives au test hypergéométrique et exporter leurs coordonnées")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dge) || is.null(opt$fai) || is.null(opt$gtf)) {
  stop("Usage: --dge DGE.csv --fai genome.fa.fai --gtf geneset.gtf [--regions regions.csv]")
}
if (!file.exists(opt$dge)) stop("DGE introuvable: ", opt$dge)
if (!file.exists(opt$fai)) stop("FAI introuvable: ", opt$fai)
if (!file.exists(opt$gtf)) stop("GTF introuvable: ", opt$gtf)
if (!is.na(opt$regions) && !file.exists(opt$regions)) stop("Régions introuvable: ", opt$regions)

# =========================
# Helpers
# =========================
out_base <- sub("\\.pdf$", "", opt$out)
writeln <- function(...) cat(paste0(..., "\n"))
write_csv <- function(df, tag) {
  fn <- sprintf("%s_%s.csv", out_base, tag)
  suppressWarnings(write.csv(df, fn, row.names = FALSE))
  writeln(sprintf("CSV écrit: %s (%d lignes)", fn, nrow(df)))
  fn
}

# Normalisation chrom
normalize_chr <- function(x) {
  x  <- as.character(x)
  xu <- toupper(x)

  out <- rep(NA_character_, length(xu))

  m <- regexpr("CHR[ _-]*(X|[1-5])", xu, perl = TRUE)
  has_chr_tag <- m > 0
  if (any(has_chr_tag)) {
    cap <- regmatches(xu, m)
    cap <- sub("^CHR[ _-]*", "", cap, perl = TRUE)
    out[has_chr_tag] <- cap
  }

  tok <- function(txt, pat) grepl(paste0("(^|[^A-Z0-9])", pat, "([^A-Z0-9]|$)"), txt, perl = TRUE)
  out[is.na(out) & tok(xu,"I")]   <- "1"
  out[is.na(out) & tok(xu,"II")]  <- "2"
  out[is.na(out) & tok(xu,"III")] <- "3"
  out[is.na(out) & tok(xu,"IV")]  <- "4"
  out[is.na(out) & tok(xu,"V")]   <- "5"
  out[is.na(out) & tok(xu,"X")]   <- "X"

  dig <- function(txt, d) grepl(paste0("(^|[^0-9])", d, "([^0-9]|$)"), txt, perl = TRUE)
  out[is.na(out) & dig(xu,"1")] <- "1"
  out[is.na(out) & dig(xu,"2")] <- "2"
  out[is.na(out) & dig(xu,"3")] <- "3"
  out[is.na(out) & dig(xu,"4")] <- "4"
  out[is.na(out) & dig(xu,"5")] <- "5"

  out[is.na(out) & grepl("MT|MITO|MITOCHON", xu)] <- "M"

  out
}

merge_overlaps <- function(tbl, chr_len) {
  if (nrow(tbl) == 0) return(tbl)
  tbl <- tbl[order(tbl$start, tbl$end), ]
  merged <- list(); cur <- tbl[1, ]
  if (nrow(tbl) > 1) {
    for (i in 2:nrow(tbl)) {
      if (tbl$start[i] <= cur$end) cur$end <- max(cur$end, tbl$end[i], na.rm = TRUE)
      else { merged[[length(merged)+1]] <- cur; cur <- tbl[i, ] }
    }
  }
  merged[[length(merged)+1]] <- cur
  res <- do.call(rbind, merged)
  res$start <- pmax(0L, as.integer(res$start))
  res$end   <- pmin(as.integer(chr_len), as.integer(res$end))
  res
}

# ============================================================
# ÉTAPE 01 — FAI
# ============================================================
writeln("=== ÉTAPE 01/06 : FAI ===")
fai_raw <- read.table(opt$fai, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(fai_raw) <- c("Chromosome", "Length", "Offset", "LineBases", "LineBytes")
fai_norm_all <- data.frame(ChrRaw = fai_raw$Chromosome,
                           Chr    = normalize_chr(fai_raw$Chromosome),
                           ChrLen = as.numeric(fai_raw$Length))
write_csv(fai_norm_all, "01_fai_normalization_debug")

fai <- subset(fai_norm_all, Chr %in% c("1","2","3","4","5","X"))
if (nrow(fai) == 0) {
  warning("FAI: aucun 1..5,X détecté. Fallback top-6.")
  top6 <- fai_norm_all[order(-fai_norm_all$ChrLen), ][1:min(6, nrow(fai_norm_all)), ]
  labels <- c("1","2","3","4","5","X")[seq_len(nrow(top6))]
  fai <- data.frame(Chr = labels, ChrLen = top6$ChrLen)
  write_csv(data.frame(Label=labels, ChrRaw=top6$ChrRaw, Len=top6$ChrLen), "01_fai_fallback_map")
} else {
  fai <- aggregate(ChrLen ~ Chr, data = fai, FUN = max)
  fai$Chr <- factor(fai$Chr, levels = c("1","2","3","4","5","X"))
  fai <- fai[order(fai$Chr), ]; fai$Chr <- as.character(fai$Chr)
}
write_csv(fai, "01_fai_norm")
writeln("Chromosomes retenus: ", paste(fai$Chr, collapse = ", "))

# ============================================================
# ÉTAPE 02 — Régions (input utilisateur)
# ============================================================
writeln("=== ÉTAPE 02/06 : Régions ===")
regions_highlight <- data.frame(Chr=character(), start=integer(), end=integer())
if (!is.na(opt$regions)) {
  rgn <- read.csv(opt$regions, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("Chr","startMb","endMb"))
  rgn$Chr <- normalize_chr(rgn$Chr)
  rgn$startMb <- as.numeric(rgn$startMb); rgn$endMb <- as.numeric(rgn$endMb)
  swap <- rgn$startMb > rgn$endMb
  if (any(swap)) {
    tmp <- rgn$startMb[swap]
    rgn$startMb[swap] <- rgn$endMb[swap]
    rgn$endMb[swap] <- tmp
  }
  rgn$start <- as.integer(round(rgn$startMb * 1e6))
  rgn$end   <- as.integer(round(rgn$endMb   * 1e6))
  rgn <- rgn[, c("Chr","start","end")]
  rgn <- merge(rgn, fai, by = "Chr", all.x = TRUE)
  if (any(is.na(rgn$ChrLen))) stop("Chromosomes inconnus dans --regions: ", paste(unique(rgn$Chr[is.na(rgn$ChrLen)]), collapse=", "))
  rgn$start <- pmax(0L, rgn$start); rgn$end <- pmin(rgn$ChrLen, rgn$end)
  rgn <- rgn[rgn$end > rgn$start, c("Chr","start","end")]

  regions_highlight <- rgn[order(rgn$Chr, rgn$start), ]
  rownames(regions_highlight) <- NULL

  write_csv(regions_highlight, "02_regions_merged_bp")
  bed_out <- sprintf("%s_02_regions.bed", out_base)
  if (nrow(regions_highlight)>0) {
    write.table(regions_highlight[, c("Chr","start","end")], bed_out,
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    writeln("BED écrit: ", bed_out)
  }
} else {
  writeln("Pas de fichier --regions fourni.")
}

# ============================================================
# ÉTAPE 03 — GTF
# ============================================================
writeln("=== ÉTAPE 03/06 : GTF ===")
gtf_data <- read.table(opt$gtf, header = FALSE, sep = "\t",
                       comment.char = "#", stringsAsFactors = FALSE)
gene_data <- gtf_data[gtf_data$V3 == "gene", ]
gene_data$GeneID <- substring(gene_data$V9, 9, 21)  # adapté à ton GTF
gene_data$Midpoint <- gene_data$V4 + (gene_data$V5 - gene_data$V4) / 2
gene_data$V1 <- normalize_chr(gene_data$V1)
gene_data <- gene_data[gene_data$V1 %in% c("1","2","3","4","5","X"), ]
gtf_prepped <- data.frame(Chr = gene_data$V1,
                          GenId = gene_data$GeneID,
                          GenMidpoint = gene_data$Midpoint)
gtf_prepped <- gtf_prepped[order(gtf_prepped$Chr,
                                 gtf_prepped$GenMidpoint,
                                 gtf_prepped$GenId), ]
write_csv(gtf_prepped, "03_gtf_genes_midpoints")
gtf_chr_stats <- aggregate(GenId ~ Chr, data = gtf_prepped, FUN = length)
names(gtf_chr_stats)[2] <- "n_genes"
write_csv(gtf_chr_stats, "03_genes_per_chr")

# ============================================================
# ÉTAPE 04 — DGE
# ============================================================
writeln("=== ÉTAPE 04/06 : DGE ===")
dge_raw <- read.csv(opt$dge, header = TRUE, check.names = FALSE)
colnames(dge_raw)[1] <- "GeneID"
names(dge_raw)[names(dge_raw) == "GeneID"] <- "GenId"
names(dge_raw)[names(dge_raw) == "log2FoldChange"] <- "log2FC"
if (!all(c("GenId","log2FC","padj") %in% names(dge_raw))) {
  stop("DGE doit contenir: GenId, log2FC (ou log2FoldChange), padj")
}
dge_clean <- dge_raw[, c("GenId","log2FC","padj")]
dge_clean <- dge_clean[order(dge_clean$GenId), ]
write_csv(dge_clean, "04_dge_clean")

# ============================================================
# ÉTAPE 05 — Merge & Filtrage
# ============================================================
writeln("=== ÉTAPE 05/06 : Merge & Filtrage ===")
fai_use <- read.csv(sprintf("%s_01_fai_norm.csv", out_base), stringsAsFactors = FALSE)
gtf_use <- read.csv(sprintf("%s_03_gtf_genes_midpoints.csv", out_base), stringsAsFactors = FALSE)
dge_use <- read.csv(sprintf("%s_04_dge_clean.csv", out_base), stringsAsFactors = FALSE)

df <- merge(gtf_use, dge_use, by = "GenId", all.x = TRUE)
df <- df[!is.na(df$padj), ]
df$padj[df$padj == 0] <- 1e-10
df$log10padj <- -log10(df$padj)
df$keep <- with(df, padj < opt$padj & (log2FC >= opt$lfc | log2FC <= -opt$lfc))
write_csv(df, "05_merged_all_genes")
write_csv(subset(df, keep), "05_deg_only")
deg_chr <- aggregate(keep ~ Chr, data = df, FUN = sum)
colnames(deg_chr)[2] <- "n_DEG"
write_csv(deg_chr, "05_DEG_per_chr")
writeln(sprintf("Total gènes filtrables: %d | DEG: %d", nrow(df), sum(df$keep)))

# ============================================================
# ÉTAPE 06 — Hypergéométrique (+ régions significatives)
# ============================================================
writeln("=== ÉTAPE 06/06 : Hypergéométrique ===")
hyper_all <- list()
for (chr in fai_use$Chr) {
  dfc <- subset(df, Chr == chr)
  end <- fai_use$ChrLen[fai_use$Chr == chr][1]
  taille_fenetre <- as.integer(opt$win)
  pas <- as.integer(opt$step)
  if (end < taille_fenetre) next

  bornes_debut <- seq(0, max(0, end - taille_fenetre), by = pas)
  bornes_fin <- bornes_debut + taille_fenetre
  fen <- data.frame(debut = bornes_debut, fin = bornes_fin)
  fen$midwindow <- (fen$fin + fen$debut) / 2
  fen$nbtot <- NA_integer_
  fen$nbdeg <- NA_integer_

  for (i in seq_len(nrow(fen))) {
    dfctot <- dfc[dfc$GenMidpoint >= fen$debut[i] &
                  dfc$GenMidpoint <= fen$fin[i], ]
    dfcdeg <- dfctot[dfctot$keep == TRUE, ]
    fen$nbtot[i] <- nrow(dfctot)
    fen$nbdeg[i] <- nrow(dfcdeg)
  }

  N_chr <- nrow(dfc)
  K_chr <- sum(dfc$keep)
  fen$HypergeomPval <- 1
  ok <- (fen$nbtot > 0) & (N_chr > 0)
  x <- fen$nbdeg[ok]
  m <- K_chr
  n <- N_chr - K_chr
  k <- fen$nbtot[ok]
  fen$HypergeomPval[ok] <- phyper(x - 1L, m, n, k,
                                  lower.tail = FALSE, log.p = FALSE)
  res <- fen
  res$Chr <- chr
  res$HypergeomFDR <- p.adjust(res$HypergeomPval, method = "BH")
  res$log10HypergeomFDR <- -log10(res$HypergeomFDR)
  hyper_all[[length(hyper_all)+1]] <- res
}

hyper_sig_regions <- data.frame(Chr=character(), start=integer(), end=integer())

if (length(hyper_all) > 0) {
  hyper_df <- do.call(rbind, hyper_all)
  hyper_df <- hyper_df[order(hyper_df$Chr, hyper_df$midwindow), ]
  write_csv(hyper_df, "06_hypergeom_windows")

  # 1) runs de fenêtres consécutives FDR < padj -> régions brut (debut[s]..fin[e])
  run_regions <- list()

  for (chr in unique(hyper_df$Chr)) {
    h_chr <- hyper_df[hyper_df$Chr == chr, ]
    is_sig <- is.finite(h_chr$HypergeomFDR) &
              !is.na(h_chr$HypergeomFDR) &
              h_chr$HypergeomFDR < opt$padj

    if (!any(is_sig)) next

    r <- rle(is_sig)
    ends   <- cumsum(r$lengths)
    starts <- c(1, head(ends, -1) + 1)

    for (k in seq_along(r$values)) {
      if (!r$values[k]) next
      s <- starts[k]
      e <- ends[k]

      region_start <- h_chr$debut[s]
      region_end   <- h_chr$fin[e]

      run_regions[[length(run_regions) + 1]] <- data.frame(
        Chr   = chr,
        start = region_start,
        end   = region_end
      )
    }
  }

  # 2) fusion des régions qui se chevauchent, par chromosome
  if (length(run_regions) > 0) {
    reg_all <- do.call(rbind, run_regions)
    reg_all <- reg_all[order(reg_all$Chr, reg_all$start), ]

    merged_all <- list()
    for (chr in unique(reg_all$Chr)) {
      chr_len <- fai_use$ChrLen[fai_use$Chr == chr][1]
      tmp <- reg_all[reg_all$Chr == chr, ]
      m   <- merge_overlaps(
        data.frame(start = as.integer(tmp$start),
                   end   = as.integer(tmp$end)),
        chr_len = chr_len
      )
      if (nrow(m) > 0) {
        m$Chr <- chr
        m <- m[, c("Chr","start","end")]
        merged_all[[length(merged_all)+1]] <- m
      }
    }

    if (length(merged_all) > 0) {
      hyper_sig_regions <- do.call(rbind, merged_all)
      hyper_sig_regions <- hyper_sig_regions[order(hyper_sig_regions$Chr,
                                                   hyper_sig_regions$start), ]
      rownames(hyper_sig_regions) <- NULL
    }
  }

} else {
  write_csv(data.frame(), "06_hypergeom_windows")
  writeln("Hypergéométrique: aucune fenêtre calculée.")
}

# Export des régions significatives (CSV + BED) – coordonnées = fenêtres fusionnées
write_csv(hyper_sig_regions, "06_hypergeom_sig_regions_merged_bp")
if (nrow(hyper_sig_regions) > 0) {
  bed_sig <- sprintf("%s_06_hypergeom_sig_regions.bed", out_base)
  write.table(hyper_sig_regions[, c("Chr","start","end")], bed_sig,
              quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  writeln("BED hypergéom significatif écrit: ", bed_sig)
}

# ============================================================
# FIGURE — Facets (1 plot, 2 tracks, colonnes par chromosome)
# ============================================================
writeln("\n=== FIGURE (facet_grid) ===")

fai_plot   <- read.csv(sprintf("%s_01_fai_norm.csv", out_base), stringsAsFactors = FALSE)
df_plot    <- read.csv(sprintf("%s_05_merged_all_genes.csv", out_base), stringsAsFactors = FALSE)
hyper_plot <- read.csv(sprintf("%s_06_hypergeom_windows.csv", out_base), stringsAsFactors = FALSE)

regions_plot <- data.frame(Chr=character(), start=integer(), end=integer())
if (!is.na(opt$regions)) {
  rp <- sprintf("%s_02_regions_merged_bp.csv", out_base)
  if (file.exists(rp)) regions_plot <- read.csv(rp, stringsAsFactors = FALSE)
}

hyper_sig_plot <- data.frame(Chr=character(), start=integer(), end=integer())
hsp <- sprintf("%s_06_hypergeom_sig_regions_merged_bp.csv", out_base)
if (file.exists(hsp)) {
  tmp_hsp <- read.csv(hsp, stringsAsFactors = FALSE)
  if (nrow(tmp_hsp) > 0) hyper_sig_plot <- tmp_hsp
}

chr_levels <- c("1","2","3","4","5","X")
chr_labs   <- setNames(as.character(1:20), as.character(1:20))
chr_labs["X"] <- "X"

df_plot$Chr    <- factor(df_plot$Chr, levels = chr_levels)
hyper_plot$Chr <- if (nrow(hyper_plot)) factor(hyper_plot$Chr, levels = chr_levels) else hyper_plot$Chr
regions_plot$Chr <- if (nrow(regions_plot)) factor(regions_plot$Chr, levels = chr_levels) else regions_plot$Chr
hyper_sig_plot$Chr <- if (nrow(hyper_sig_plot)) factor(hyper_sig_plot$Chr, levels = chr_levels) else hyper_sig_plot$Chr
fai_plot$Chr <- factor(fai_plot$Chr, levels = chr_levels)

minlog10padj <- if (sum(df_plot$keep)>0) {
  min(df_plot$log10padj[df_plot$keep], na.rm=TRUE)
} else {
  min(df_plot$log10padj, na.rm=TRUE)
}
minlog10padj <- ifelse(is.finite(minlog10padj), minlog10padj, 1)

genes_df <- df_plot %>%
  dplyr::filter(keep,
                is.finite(log2FC),
                is.finite(log10padj),
                is.finite(GenMidpoint)) %>%
  dplyr::mutate(
    log2FC   = pmax(pmin(log2FC,  4), -4),
    size_var = pmin(pmax(log10padj, minlog10padj), 30),
    track    = "DEG"
  )

enrich_df <- if (nrow(hyper_plot)) {
  hyper_plot %>%
    dplyr::mutate(track = "Enrichment") %>%
    dplyr::filter(is.finite(log10HypergeomFDR))
} else hyper_plot

regions_df <- if (nrow(regions_plot)) {
  rbind(
    transform(regions_plot, track = "DEG"),
    transform(regions_plot, track = "Enrichment")
  )
} else regions_plot

regions_hyper_df <- if (nrow(hyper_sig_plot) && isTRUE(opt$hyperHighlight)) {
  rbind(
    transform(hyper_sig_plot, track = "DEG"),
    transform(hyper_sig_plot, track = "Enrichment")
  )
} else hyper_sig_plot[0, ]

if (nrow(enrich_df) > 0) {
  enr_M <- max(enrich_df$log10HypergeomFDR, na.rm = TRUE)
  enr_ymax <- if (is.finite(enr_M)) enr_M + 1 else 1
} else enr_ymax <- 1
enr_ymax <- ceiling(enr_ymax)

deg_breaks_fixed  <- seq(-4, 4, by = 2)
enr_breaks_every2 <- seq(0, enr_ymax, by = 2)

lab_chr <- ggplot2::as_labeller(chr_labs)
lab_mb  <- scales::label_number(accuracy = 5, scale = 1e-6)

lab_track <- as_labeller(
  c(DEG = "-log[2](FC)", Enrichment = "-log[10](FDR)"),
  default = label_parsed
)

hline_deg  <- data.frame(track = "DEG",        y = 0)
hline_enr  <- data.frame(track = "Enrichment", y = -log10(opt$padj))

p <- ggplot() +
  { if (nrow(regions_df) > 0)
      geom_rect(
        data = regions_df,
        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "indianred1",
        alpha = 0.25
      ) else NULL } +
  { if (nrow(regions_hyper_df) > 0 && isTRUE(opt$hyperHighlight))
      geom_rect(
        data = regions_hyper_df,
        aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "indianred1",
        alpha = 0.45
      ) else NULL } +
  { if (nrow(genes_df) > 0)
      geom_point(
        data = genes_df,
        aes(x = GenMidpoint, y = log2FC, size = size_var, color = log2FC >= 0),
        stroke = 0, alpha = 0.8, na.rm = TRUE
      ) else NULL } +
  { if (nrow(enrich_df) > 0)
      geom_line(
        data = enrich_df,
        aes(x = midwindow, y = log10HypergeomFDR, group = 1),
        linewidth = 0.4, color = "blue"
      ) else NULL } +
  geom_hline(data = hline_deg,
             aes(yintercept = y),
             linewidth = 0.6,
             linetype = "solid",
             color = "grey40") +
  geom_hline(data = hline_enr,
             aes(yintercept = y),
             linewidth = 0.6,
             linetype = "dashed",
             color = "red") +
  scale_color_manual(values = c("TRUE"="orange","FALSE"="blue"), guide = "none") +
  scale_size_continuous(
    name = expression("-log"[10]*"p"[adj]),
    limits = c(minlog10padj, 30),
    breaks = c(minlog10padj, 10, 30),
    labels = round(c(minlog10padj, 10, 30), 1),
    range  = c(1.2, 4.5)
  ) +
  facet_grid(
    rows     = vars(track),
    cols     = vars(Chr),
    scales   = "free",
    space    = "free_x",
    labeller = labeller(Chr = lab_chr, track = lab_track),
    switch   = "y"
  ) +
  scale_x_continuous(labels = lab_mb) +
  coord_cartesian(clip = "off") +
  labs(title = opt$title, x = "Genomic position (Mb)", y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 20),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = "right",
    strip.placement = "outside",
    strip.background= element_blank(),
    strip.text.x    = element_text(size = 18),
    strip.text.y    = element_text(size = 14),
    axis.title.x    = element_text(size = 18),
    axis.text.x     = element_text(size = 16),
    axis.text.y     = element_text(size = 16),
    legend.title    = element_text(size = 16),
    legend.text     = element_text(size = 15),
    panel.grid.major= element_line(colour = "grey85", linewidth = 0.5),
    panel.grid.minor= element_blank()
  )

if (requireNamespace("ggh4x", quietly = TRUE)) {
  p <- p +
    ggh4x::facetted_pos_scales(
      y = list(
        "DEG" = ggplot2::scale_y_continuous(
          limits = c(-4.5, 4.5),
          breaks = deg_breaks_fixed
        ),
        "Enrichment" = ggplot2::scale_y_continuous(
          limits = c(0, enr_ymax),
          breaks = enr_breaks_every2
        )
      )
    )
} else {
  warning("ggh4x non disponible : Y global appliqué.")
  p <- p + scale_y_continuous(
    limits = c(-4.5, max(4.5, enr_ymax)),
    breaks = seq(-4, max(4, enr_ymax), by = 2)
  )
}

g <- ggplotGrob(p)
panel_rows <- unique(g$layout$t[grepl("^panel", g$layout$name)])
if (length(panel_rows) >= 2) {
  g$heights[panel_rows[1]] <- unit(3, "null")
  g$heights[panel_rows[2]] <- unit(2, "null")
}

ggsave(opt$out, plot = g, width = opt$width, height = opt$height)
writeln(sprintf("PDF écrit: %s", opt$out))
invisible(NULL)
