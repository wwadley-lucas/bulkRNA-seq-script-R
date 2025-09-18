########################
## GSEA FROM SCRATCH  ##
## MANUAL FILE INPUTS ##
########################

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(fgsea)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
  library(ragg)
  library(cowplot)
})

# ==============================
# 1) MANUAL PATHS + PARAMETERS
# ==============================

# ---- Input files ----
expr_file <- "PATH/TO/FILE/EXPRESSION FILE"
cls_file  <- "PATH/TO/FILE/GSEA CLS FILE"

# ---- Groups to compare (must match CLS labels) ---- 
#manual inputs
groupA <- "GROUP 1"
groupB <- "GROUP 2"

# ---- Species and MSigDB collection ----
species         <- "Mus musculus"   # or "Homo sapiens"
msig_collection <- "C2"              # e.g. "H", "C2", "C5", ...
msig_subcol     <- "CP:BIOCARTA"             # e.g. "GO:BP" when collection="C5", else NULL

# ---- Output and plotting ----
out_dir    <- "PATH/TO/FILE/OUTPUT DIRECTORY"
top_plots  <- 100                   # hard cap on plots to make
pval_cut   <- 0.05                  # plot sets with raw p-value <= this
NES_positive_group <- "A"           # "A" => NES>0 enriched in A; set "B" to flip

# ==============================
# 2) HELPERS
# ==============================

read_cls <- function(path) {
  lines <- readLines(path, warn = FALSE)
  stopifnot(length(lines) >= 3)
  classes <- strsplit(trimws(sub("^#\\s*", "", lines[2])), "\\s+")[[1]]
  labels  <- strsplit(trimws(lines[3]), "\\s+")[[1]]
  list(classes = classes[nzchar(classes)], labels = labels[nzchar(labels)])
}

read_matrix <- function(path) {
  dt <- fread(path)
  genes <- dt[[1]]
  mat <- as.data.frame(dt[, -1])
  rownames(mat) <- genes
  mat
}

# Broad-style Signal2Noise (USE_MEDIAN = FALSE, FIX_LOW = TRUE)
# s2n = (meanA - meanB) / (sdA + sdB), with sd floored at 0.2; no epsilon
signal_to_noise <- function(X, labs, A, B, sd_floor = 0.2) {
  ia <- which(labs == A); ib <- which(labs == B)
  stopifnot(length(ia) >= 2, length(ib) >= 2)
  
  mA <- rowMeans(X[, ia, drop=FALSE])
  mB <- rowMeans(X[, ib, drop=FALSE])
  sA <- pmax(apply(X[, ia, drop=FALSE], 1, sd), sd_floor)
  sB <- pmax(apply(X[, ib, drop=FALSE], 1, sd), sd_floor)
  
  s2n <- (mA - mB) / (sA + sB)
  names(s2n) <- rownames(X)
  s2n <- s2n[is.finite(s2n)]
  
  # collapse duplicate gene IDs by max |score|; uppercase to match MSigDB
  ranked <- vapply(split(s2n, toupper(names(s2n))),
                   function(v) v[which.max(abs(v))], numeric(1))
  ranked <- sort(ranked, decreasing = TRUE)
  stopifnot(is.numeric(ranked), length(ranked) > 50L, !is.null(names(ranked)),
            all(diff(unname(ranked)) <= 0))
  ranked
}

# TERM2GENE from MSigDB using collection + optional subcollection
msig_TERM2GENE <- function(species, collection, subcollection = NULL) {
  msig <- msigdbr(
    species       = species,
    collection    = collection,
    subcollection = subcollection
  )
  df <- unique(msig[, c("gs_name", "gene_symbol")])
  df$gs_name      <- as.character(df$gs_name)
  df$gene_symbol  <- toupper(as.character(df$gene_symbol))  # match our ranks
  df
}

save_png_white <- function(filename, plot, width = 7.6, height = 5.8, dpi = 300) {
  ragg::agg_png(filename, width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot); dev.off()
}

# Make ES / hits / ranked panels align & look consistent
normalize_panel <- function(g, add_border = TRUE) {
  g +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      plot.margin  = margin(4, 6, 2, 6),
      panel.border = if (add_border) element_rect(fill = NA, colour = "black", linewidth = 0.5)
      else element_blank()
    )
}

# ==============================
# 3) MAIN
# ==============================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(out_dir, sprintf("%s_vs_%s_plots", groupA, groupB))
dir.create(plot_dir, showWarnings = FALSE)

message("Reading expression: ", expr_file)
# After reading the raw counts matrix X
X <- read_matrix(expr_file)
# Transform raw counts -> log2(CPM + 1)
library(edgeR)
X <- cpm(X, log = TRUE, prior.count = 1)  # or log2((counts/sum)*1e6 + 1)

message("Reading class file: ", cls_file)
cls <- read_cls(cls_file)

# CLS labels must equal sample columns
nsamp <- ncol(X)
if (length(cls$labels) != nsamp) {
  stop(sprintf("CLS labels (%d) != # samples in matrix (%d).\nCheck expression columns (excluding gene column) vs CLS labels.",
               length(cls$labels), nsamp))
}

# Ensure requested groups exist
if (!all(c(groupA, groupB) %in% cls$labels)) {
  stop(sprintf("groupA='%s' or groupB='%s' not found in CLS labels. Available: %s",
               groupA, groupB, paste(sort(unique(cls$labels)), collapse = ", ")))
}

# (Optional) Transform: Broad GSEA typically uses already log-normalized inputs.
# If your matrix is raw counts, normalize/log outside of this script.
X <- as.matrix(X)

# Rank genes (POS => enriched in A)
stats <- signal_to_noise(X, cls$labels, groupA, groupB)

# Optional: flip so NES>0 == enriched in B (and keep decreasing order)
if (identical(NES_positive_group, "B")) {
  stats <- sort(-stats, decreasing = TRUE)
}

# TERM2GENE from chosen collection/subcollection
t2g      <- msig_TERM2GENE(species, msig_collection, msig_subcol)
pathways <- split(t2g$gene_symbol, t2g$gs_name)

# Overlap sanity
overlaps <- vapply(pathways, function(g) sum(g %in% names(stats)), integer(1))
message("Pathways with >=10 overlapping genes: ",
        sum(overlaps >= 10), " / ", length(pathways))
if (!any(overlaps >= 10)) {
  stop("No gene sets overlap your ranked list. Check gene IDs (SYMBOLs) and species/collection.")
}

# ---------- fgsea table (settings to mirror Broad) ----------
set.seed(as.integer(Sys.time()))
fg <- fgsea::fgseaMultilevel(
  pathways    = pathways,
  stats       = stats,
  minSize     = 5,       # was 10
  maxSize     = 1000,    # was 5000
  gseaParam   = 0,       # classic (unweighted) instead of 1
  nPermSimple = 1000,
  eps         = 0
)
fg <- fg[order(fg$pval, -fg$NES), ]
fg <- fg[!is.na(fg$NES), ]  # drop any sets with NA NES before plotting

out_csv <- file.path(out_dir, sprintf("%s_vs_%s__fgsea_results.csv", groupA, groupB))
fwrite(as.data.table(fg), out_csv)
message("Saved fgsea table: ", out_csv)

# ---------- clusterProfiler object (for pretty ES curve) ----------
gsea_cp <- clusterProfiler::GSEA(
  geneList      = stats,
  TERM2GENE     = t2g,
  exponent      = 0,     # classic (unweighted) instead of default 1
  minGSSize     = 5,
  maxGSSize     = 1000,
  pvalueCutoff  = 1.0,
  pAdjustMethod = "BH",
  seed          = TRUE,
  verbose       = FALSE
)

# ---------- choose sets to plot (by raw p-value) ----------
pick <- subset(fg, !is.na(pval) & pval <= pval_cut)
if (nrow(pick) == 0) {
  message("No pathways pass pval <= ", pval_cut, "; no plots will be generated.")
} else {
  pick <- pick[order(pick$pval, -pick$NES), , drop = FALSE]
  if (nrow(pick) > top_plots) pick <- head(pick, top_plots)
  message("Plotting ", nrow(pick), " pathways with pval <= ", pval_cut)
}

# ---------- plotting loop (aligned panels + dashed zero lines) ----------
for (i in seq_len(nrow(pick))) {
  gs <- pick$pathway[i]
  
  title_main <- stringr::str_wrap(
    paste0(gsub("^HALLMARK_|^GOBP_|^GO_BP_", "", gs),
           "  (NES=", sprintf("%.2f", pick$NES[i]),
           ", P=", format(pick$pval[i], digits = 2), ")"),
    width = 60
  )
  title_hint <- paste0("Left: Enriched in ", groupA, "   |   Right: Enriched in ", groupB)
  full_title <- paste(title_main, title_hint, sep = "\n")
  
  # ES (subplot 1) with dashed 0-line
  p_es <- enrichplot::gseaplot2(
    gsea_cp, geneSetID = gs, base_size = 11, subplots = 1
  ) +
    ggplot2::labs(title = full_title) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    ggplot2::theme(plot.title = element_text(hjust = 0, face = "bold"))
  p_es <- normalize_panel(p_es)
  
  # Hit indices (subplot 2)
  p_hits <- enrichplot::gseaplot2(
    gsea_cp, geneSetID = gs, base_size = 11, subplots = 2
  )
  p_hits <- normalize_panel(p_hits)
  
  # Ranked metric (subplot 3) with dashed 0-line
  p_rank <- enrichplot::gseaplot2(
    gsea_cp, geneSetID = gs, base_size = 11, subplots = 3
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed")
  p_rank <- normalize_panel(p_rank)
  
  # Align and stack
  aligned <- cowplot::align_plots(p_es, p_hits, p_rank, align = "v", axis = "lr")
  p_final <- cowplot::plot_grid(plotlist = aligned, ncol = 1,
                                rel_heights = c(1.1, 0.30, 1.0))
  
  out_png <- file.path(
    plot_dir,
    paste0(groupA, "_vs_", groupB, "__", gsub("[^A-Za-z0-9]+","_", gs), ".png")
  )
  save_png_white(out_png, p_final, width = 7.6, height = 5.8, dpi = 300)
}

message("Plots saved to: ", plot_dir)
message("DONE: ", groupA, " vs ", groupB)
message(sprintf("MSigDB: collection=%s  subcollection=%s  species=%s",
                msig_collection, ifelse(is.null(msig_subcol), "NULL", msig_subcol), species))
message(sprintf("Contrast: %s vs %s  (NES>0 → %s)", groupA, groupB, NES_positive_group))