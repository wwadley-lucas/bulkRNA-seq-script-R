##### PACKAGES #####
suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(GSVA)
  library(fgsea)
  library(ggpubr)
  library(stringr)
  library(org.Mm.eg.db)  
  library(org.Hs.eg.db)# change to org.Hs.eg.db if human
  library(AnnotationDbi)
  library(clusterProfiler)
  library(enrichplot)
  library(cowplot)
  library(stringr)
  library(patchwork)
  library(grid)
  library(ragg)
})

##### CONFIG: designs (edit here) #####
design_formula_default <- ~ Exposure
design_formula_grouped <- ~ Exposure

##### CONFIG: paths & options #####
paths <- list(
  project_dir  = "PATH/TO/FILE/DIRECTORY",
  counts_csv   = "Tables/All probes.csv", #expression raw counts
  meta_csv     = "Tables/metaData.csv",
  meta_grouped = "Tables/metaData.csv",  # reuse if only one meta
  out_base     = "Figures/"
)
FIG <- list(W = 8, H = 6, DPI = 600)
species_msigdb <- "Mus musculus"   # or "Homo sapiens"
##### THRESHOLDS #####
thr <- list(
  pval          = 0.05,
  log2fc        = 1,      # for quick DEG heatmaps
  volcano_p     = 0.05,
  volcano_log2  = 1,
  go_alpha_padj = 0.05,
  go_lfc_cut    = 2
)

##### VOLCANO SETTINGS #####
volcano <- list(group_col = "Group", y_col = "pvalue", W = 7, H = 6, DPI = 600)

##### GO SETTINGS #####
GO <- list(
  group_col     = "Group",
  minGSSize     = 15,
  maxGSSize     = 500,
  do_simplify   = TRUE,
  min_genes_dir = 5,
  ont_choice    = "ALL",   # "BP"/"MF"/"CC"/"ALL"
  W = 9, H = 7.5, DPI = 300
)
##### OUTPUT SUBFOLDERS #####
# Define all output subfolders
outs <- list(
  plots               = file.path(paths$out_base, "plots"),
  tables              = file.path(paths$out_base, "tables"),
  heatmaps_all        = file.path(paths$out_base, "heatmaps_all"),
  heatmaps_sig        = file.path(paths$out_base, "heatmaps_sig"),
  volcano_plots       = file.path(paths$out_base, "volcanoes/plots"),
  volcano_tables      = file.path(paths$out_base, "volcanoes/tables"),
  go_tables_raw       = file.path(paths$out_base, "GO/tables/raw"),
  go_tables_simpl     = file.path(paths$out_base, "GO/tables/simpl"),
  go_tables_cc        = file.path(paths$out_base, "GO/tables/compareCluster"),
  go_plots_bar_raw    = file.path(paths$out_base, "GO/plots/bar_raw"),
  go_plots_bar_simpl  = file.path(paths$out_base, "GO/plots/bar_simpl"),
  go_plots_dot_raw    = file.path(paths$out_base, "GO/plots/dot_raw"),
  go_plots_dot_simpl  = file.path(paths$out_base, "GO/plots/dot_simpl"),
  go_plots_cc         = file.path(paths$out_base, "GO/plots/compareCluster"),
  gsea_tables         = file.path(paths$out_base, "GSEA/tables"),
  gsea_plots          = file.path(paths$out_base, "GSEA/plots")
)

# Create all directories in one go
dir_create <- function(dirs) {
  invisible(lapply(dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)))
}

dir_create(outs)  # <- pass the whole list, not do.call

##### HELPERS #####
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)
wrap_lab  <- function(x, width = 50) stringr::str_wrap(x, width)

save_plot_pdf <- function(path, plot, w = FIG$W, h = FIG$H) {
  ggplot2::ggsave(path, plot = plot, width = w, height = h, device = cairo_pdf)
}
save_plot_png <- function(path, plot, w = FIG$W, h = FIG$H, dpi = FIG$DPI) {
  ggplot2::ggsave(path, plot = plot, width = w, height = h, dpi = dpi, device = ragg::agg_png)
}
save_plot_both <- function(stem, plot, w = FIG$W, h = FIG$H, dpi = FIG$DPI) {
  save_plot_pdf(paste0(stem, ".pdf"), plot, w, h)
  save_plot_png(paste0(stem, ".png"), plot, w, h, dpi)
}

wrap_go_labels <- function(p, width = 40, y_text_size = 6, left_margin = 80) {
  p +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = width)) +
    theme(axis.text.y = element_text(size = y_text_size),
          plot.margin = margin(8, 14, 8, left_margin))
}
shrink_dots <- function(p, range = c(1.2, 4), title = "Count") {
  p + scale_size_continuous(range = range) + guides(size = guide_legend(title = title))
}
map_to_symbol <- function(ids) {
  if (length(ids) == 0) return(character())
  ids_clean <- sub("\\..*$", "", ids)
  ann <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(ids_clean),
                               keytype = "ENSEMBL", columns = c("SYMBOL"))
  ann <- ann[!is.na(ann$SYMBOL) & !duplicated(ann$ENSEMBL), ]
  symbols <- ann$SYMBOL[match(ids_clean, ann$ENSEMBL)]
  names(symbols) <- ids
  symbols
}
ensure_symbols <- function(id_vec) {
  if (length(id_vec) == 0) return(character())
  looks_ens <- any(grepl("^ENSMUSG", id_vec))
  if (!looks_ens) return(id_vec)
  sym <- map_to_symbol(id_vec); out <- sym[!is.na(sym)]; unname(out)
}

# design utils
design_vars <- function(fml) unique(all.vars(fml))
check_vars_in_coldata <- function(cd, vars) {
  ok <- vars[vars %in% colnames(cd)]
  if (length(ok) < length(vars)) warning("Missing in colData: ", paste(setdiff(vars, ok), collapse = ", "))
  ok
}

# PCA that adapts color/shape/fill + ellipse rules you asked for
plot_adaptive_pca <- function(dds, design_fml, out_stem, ntop = 2000, title = "PCA") {
  vsd  <- vst(dds, blind = FALSE)
  vars <- unique(all.vars(design_fml))
  vars <- vars[vars %in% colnames(colData(dds))]
  if (length(vars) == 0) stop("No design variables from the formula are present in colData.")
  for (v in vars) if (!is.factor(colData(vsd)[[v]])) colData(vsd)[[v]] <- factor(colData(vsd)[[v]])
  use_vars <- vars[seq_len(min(3, length(vars)))]
  pdat     <- plotPCA(vsd, intgroup = use_vars, ntop = ntop, returnData = TRUE)
  pv       <- round(100 * attr(pdat, "percentVar"))
  k <- length(use_vars)
  if (k == 1) {
    v1 <- use_vars[1]
    p <- ggplot(pdat, aes(PC1, PC2, color = .data[[v1]])) +
      geom_point(size = 3) +
      ggforce::geom_mark_ellipse(aes(fill = .data[[v1]], group = .data[[v1]]),
                                 alpha = .20, expand = unit(10, "mm"), show.legend = FALSE)
  } else if (k == 2) {
    v1 <- use_vars[1]; v2 <- use_vars[2]
    p <- ggplot(pdat, aes(PC1, PC2, color = .data[[v1]], shape = .data[[v2]])) +
      geom_point(size = 3) +
      ggforce::geom_mark_ellipse(aes(group = interaction(.data[[v1]], .data[[v2]])),
                                 alpha = .08, fill = NA, expand = unit(10, "mm"), show.legend = FALSE)
  } else {
    v1 <- use_vars[1]; v2 <- use_vars[2]; v3 <- use_vars[3]
    p <- ggplot(pdat, aes(PC1, PC2, color = .data[[v1]], shape = .data[[v2]])) +
      geom_point(size = 3) +
      ggforce::geom_mark_ellipse(aes(fill = .data[[v3]],
                                     group = interaction(.data[[v1]], .data[[v2]], .data[[v3]])),
                                 alpha = .20, expand = unit(10, "mm"), show.legend = FALSE)
  }
  p <- p +
    xlab(paste0("PC1: ", pv[1], "% variance")) +
    ylab(paste0("PC2: ", pv[2], "% variance")) +
    ggtitle(title) +
    coord_fixed() +
    theme_bw() +
    theme(panel.border = element_rect(linewidth = 2))
  save_plot_both(out_stem, p, w = 8, h = 6, dpi = 300)
  invisible(p)
}

# Title helper for GO plots
add_title <- function(p, title, subtitle = NULL, width = 60) {
  p +
    labs(title = stringr::str_wrap(title, width),
         subtitle = if (!is.null(subtitle)) stringr::str_wrap(subtitle, width) else NULL) +
    theme(plot.title.position = "plot",
          plot.title   = element_text(lineheight = 0.98, margin = margin(b = 6)),
          plot.subtitle= element_text(lineheight = 0.98, margin = margin(b = 8)))
}

if (nzchar(paths$project_dir)) setwd(paths$project_dir)

##### LOAD COUNTS/META #####
df <- read.csv(paths$counts_csv, header = TRUE, sep = ",", row.names = 1)
countData <- as.matrix(df)
countData[is.na(countData) | is.nan(countData) | is.infinite(countData)] <- 0
countData <- countData[apply(countData, 1, var) != 0, , drop = FALSE]

metaData         <- read.csv(paths$meta_csv, header = TRUE, sep = ",", row.names = 1, fill = TRUE)
metaData_grouped <- read.csv(paths$meta_grouped, header = TRUE, sep = ",", row.names = 1, fill = TRUE)

stopifnot(identical(colnames(countData), rownames(metaData)))
stopifnot(identical(colnames(countData), rownames(metaData_grouped)))

##### CREATE DESEQ2 OBJS #####
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = design_formula_default)
dds <- DESeq(dds); estimateSizeFactors(dds); estimateDispersions(dds); nbinomWaldTest(dds)

dds_grouped <- DESeqDataSetFromMatrix(countData = countData, colData = metaData_grouped, design = design_formula_grouped)
dds_grouped <- DESeq(dds_grouped); estimateSizeFactors(dds_grouped); estimateDispersions(dds_grouped); nbinomWaldTest(dds_grouped)

##### AUTO-DETECT MAIN GROUPING VAR #####
main_var <- all.vars(design_formula_grouped)[1]   # e.g., "Exposure"
volcano$group_col <- main_var
GO$group_col      <- main_var

##### PRELIMINARY HEATMAPS #####
res       <- results(dds)
res_group <- results(dds_grouped)
DEG_genes          <- rownames(res)[which(res$pvalue < thr$pval & abs(res$log2FoldChange) > thr$log2fc)]
DEG_genes_grouped  <- rownames(res_group)[which(res_group$pvalue < thr$pval & abs(res_group$log2FoldChange) > thr$log2fc)]
dir.create(outs$tables, showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(gene = DEG_genes),          file.path(outs$tables, "DEG_genes.csv"),          row.names = FALSE)
write.csv(data.frame(gene = DEG_genes_grouped),  file.path(outs$tables, "DEG_genes_grouped.csv"),  row.names = FALSE)

count_matrix         <- counts(dds, normalized = TRUE)
count_matrix_grouped <- counts(dds_grouped, normalized = TRUE)
DEG_count_matrix         <- count_matrix[DEG_genes, , drop = FALSE]
DEG_count_matrix_grouped <- count_matrix_grouped[DEG_genes_grouped, , drop = FALSE]

if (nrow(DEG_count_matrix) >= 2) {
  ph1 <- pheatmap(DEG_count_matrix, cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = FALSE, show_colnames = TRUE, border_color = NA,
                  fontsize = 10, cutree_rows = 10, clustering_distance_rows = "correlation",
                  scale = "row", color = colorRampPalette(c("blue","white","red"))(100), fontsize_row = 5)
  pdf(file.path(outs$plots, "DEG_count_matrix_heatmap.pdf"), width = 8, height = 6); grid::grid.draw(ph1$gtable); dev.off()
} else message("No DEG rows (ungrouped).")

if (nrow(DEG_count_matrix_grouped) >= 2) {
  ph2 <- pheatmap(DEG_count_matrix_grouped, cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = FALSE, show_colnames = TRUE, border_color = NA,
                  fontsize = 10, cutree_rows = 10, clustering_distance_rows = "correlation",
                  scale = "row", color = colorRampPalette(c("blue","white","red"))(100), fontsize_row = 5)
  pdf(file.path(outs$plots, "DEG_count_matrix_grouped_heatmap.pdf"), width = 8, height = 6); grid::grid.draw(ph2$gtable); dev.off()
} else message("No DEG rows (grouped).")

##### MSigDB HEATMAPS #####
significant_genes <- rownames(res)[which(res$pvalue < thr$pval & abs(res$log2FoldChange) > thr$log2fc)]
gene_set <- msigdbr(species = species_msigdb, collection = "C2", subcollection = "CP:BIOCARTA")
hallmark_names <- sort(unique(gene_set$gs_name))
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

plot_one_set <- function(set_name, matT, matF) {
  specific_genes <- gene_set %>% dplyr::filter(gs_name == set_name) %>% dplyr::pull(gene_symbol) %>% unique()
  mat_all <- matT[rownames(matT) %in% specific_genes, , drop = FALSE]
  if (nrow(mat_all) >= 2) {
    pheatmap(mat_all, scale = "row", clustering_distance_rows = "correlation", clustering_method = "complete",
             cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE,
             border_color = NA, color = my_palette, fontsize_row = 10, main = set_name,
             filename = file.path(outs$heatmaps_all, paste0(safe_name(set_name), ".png")))
  }
  filtered_genes <- intersect(specific_genes, significant_genes)
  mat_sig <- matF[rownames(matF) %in% filtered_genes, , drop = FALSE]
  if (nrow(mat_sig) >= 2) {
    pheatmap(mat_sig, scale = "row", clustering_distance_rows = "correlation", clustering_method = "complete",
             cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE,
             border_color = NA, color = my_palette, fontsize_row = 10,
             main = paste0(set_name, " [p<", thr$pval, "]"),
             filename = file.path(outs$heatmaps_sig, paste0(safe_name(set_name), "_sig.png")))
  }
}
vsd_T <- vst(dds, blind = TRUE);  matT <- assay(vsd_T)
vsd_F <- vst(dds, blind = FALSE); matF <- assay(vsd_F)
invisible(lapply(hallmark_names, function(nm) { message("Plotting: ", nm); plot_one_set(nm, matT, matF) }))

##### PCA (ADAPTIVE) #####
plot_adaptive_pca(dds,         design_formula_default, out_stem = file.path(outs$plots, "PCA_RNASeq"),          title = "PCA of RNA-seq (default design)")
plot_adaptive_pca(dds_grouped, design_formula_grouped, out_stem = file.path(outs$plots, "PCA_RNASeq_grouped"),  title = "PCA of RNA-seq (grouped design)")

##### OTHER DIAGNOSTICS #####
par(mfrow = c(1,1))
plotDispEsts(dds, ylim = c(1e-6, 1e2))
hist(res$pvalue, breaks = 20, col = "grey")
qs    <- c(0, quantile(res$baseMean[res$baseMean > 0], 0:7/7))
bins  <- cut(res$baseMean, qs)
levels(bins) <- paste0("~", round(.5*qs[-1] + .5*qs[-length(qs)]))
ratios <- tapply(res$pvalue, bins, function(p) mean(p < .01, na.rm = TRUE))
barplot(ratios, xlab = "mean normalized count", ylab = "ratio of small p-values")

vsd <- vst(dds)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(rownames(colData(vsd)))
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, show_colnames = TRUE, show_rownames = TRUE)

##### VOLCANO PLOTS #####
dds_grouped <- dds_grouped[rowSums(counts(dds_grouped)) > 10, ]
if (!volcano$group_col %in% colnames(colData(dds_grouped))) {
  warning("volcano$group_col = '", volcano$group_col, "' not in colData(dds_grouped). Volcano step will be skipped.")
} else {
  groups <- levels(droplevels(colData(dds_grouped)[[volcano$group_col]]))
  if (length(groups) >= 2) {
    pairs <- combn(groups, 2, simplify = FALSE)
    plot_volcano_for <- function(A, B) {
      message("Contrast: ", A, " vs ", B, "  (log2FC>0 => higher in ", A, ")")
      res <- results(dds_grouped, contrast = c(volcano$group_col, A, B), alpha = thr$volcano_p)
      df  <- as.data.frame(res) %>% filter(!is.na(.data[[volcano$y_col]]), !is.na(log2FoldChange))
      base      <- paste0(safe_name(A), "_vs_", safe_name(B))
      title_txt <- paste0(A, " vs ", B)
      
      vp <- EnhancedVolcano(
        df, lab = rownames(df), x = "log2FoldChange", y = volcano$y_col,
        title = wrap_lab(title_txt, 40), titleLabSize = 12,
        xlab = bquote(~Log[2]~ 'fold change'),
        pCutoff = thr$volcano_p, FCcutoff = thr$volcano_log2,
        pointSize = 0.75, labSize = 4.0, colAlpha = 0.75,
        legendPosition = "none", drawConnectors = F, max.overlaps = 1000,
      ) + theme(plot.title.position = "plot",
                plot.title = element_text(lineheight = 0.95, margin = margin(b = 6)),
                plot.margin = margin(10, 18, 10, 10))
      ggsave(file.path(outs$volcano_plots, paste0(base, ".png")), vp, width = volcano$W, height = volcano$H, dpi = volcano$DPI)
      ggsave(file.path(outs$volcano_plots, paste0(base, ".pdf")), vp, width = volcano$W, height = volcano$H, device = cairo_pdf)
      write.csv(df, file.path(outs$volcano_tables, paste0(base, ".csv")), row.names = TRUE)
      invisible(TRUE)
    }
    invisible(lapply(pairs, function(p) plot_volcano_for(p[1], p[2])))
    while (!is.null(dev.list())) dev.off()
  } else {
    message("Not enough groups in '", volcano$group_col, "' for volcano plots.")
  }
}

##### ssGSEA #####
vsd_ssgsea <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd_ssgsea)

msigdb_gene_sets <- msigdbr(species = species_msigdb, collection = "H")
msigdb_gene_sets_list <- split(msigdb_gene_sets$gene_symbol, msigdb_gene_sets$gs_name)

ssgsea_params <- ssgseaParam(expr_matrix, msigdb_gene_sets_list)
ssgsea_scores <- gsva(ssgsea_params)

# save the scores table too (optional but handy)
dir.create(outs$tables, showWarnings = FALSE, recursive = TRUE)
write.csv(ssgsea_scores,
          file.path(outs$tables, "ssgsea_scores.csv"),
          row.names = TRUE)

# draw once, then save PNG + PDF into outs$plots
ph_ssg <- pheatmap(
  ssgsea_scores, scale = "row",
  clustering_distance_rows = "correlation",
  clustering_method = "complete",
  cluster_cols = FALSE, show_rownames = TRUE, fontsize_row = 4, fontsize_col = 5,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)

# PDF
pdf(file.path(outs$plots, "ssGSEA_heatmap.pdf"), width = 8, height = 6)
grid::grid.draw(ph_ssg$gtable)
dev.off()

# PNG
png(file.path(outs$plots, "ssGSEA_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
grid::grid.draw(ph_ssg$gtable)
dev.off()

##### GO ANALYSIS (ALL/BP/MF/CC) #####
# parser that respects the main variable prefix in resultsNames
parse_AB <- function(coef_name, main_var) {
  core <- sub(paste0("^", main_var, "_"), "", coef_name)
  A <- sub("_vs_.*$", "", core)
  B <- sub("^.*_vs_", "", core)
  list(A = A, B = B)
}

add_title <- function(p, title, subtitle = NULL, width = 60) {
  p + labs(title = stringr::str_wrap(title, width),
           subtitle = if (!is.null(subtitle)) stringr::str_wrap(subtitle, width) else NULL) +
    theme(plot.title.position = "plot",
          plot.title = element_text(lineheight = 0.98, margin = margin(b = 6)),
          plot.subtitle = element_text(lineheight = 0.98, margin = margin(b = 8)))
}

run_go_for_coef <- function(coef_name) {
  message("\n=== GO: ", coef_name, " ===")
  res <- results(dds_grouped, name = coef_name, alpha = thr$go_alpha_padj)
  df  <- as.data.frame(res)
  
  row_ids <- rownames(df)
  gene_list_vals <- df$log2FoldChange; names(gene_list_vals) <- row_ids
  gene_list_vals <- sort(na.omit(gene_list_vals), decreasing = TRUE)
  universe_syms <- ensure_symbols(names(gene_list_vals))
  if (length(universe_syms) < 20) { warning("Small universe for ", coef_name); return(invisible(NULL)) }
  
  sig_df <- df %>% dplyr::filter(!is.na(padj), padj < thr$go_alpha_padj, !is.na(log2FoldChange))
  sig_ids <- rownames(sig_df); if (!length(sig_ids)) { warning("No sig genes: ", coef_name); return(invisible(NULL)) }
  
  sym_map <- data.frame(id = sig_ids, symbol = ensure_symbols(sig_ids), lfc = sig_df$log2FoldChange,
                        stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(symbol)) %>% dplyr::group_by(symbol) %>%
    dplyr::slice_max(order_by = abs(lfc), n = 1, with_ties = FALSE) %>% dplyr::ungroup()
  
  genes_all <- sym_map %>% dplyr::filter(abs(lfc) > thr$go_lfc_cut) %>% dplyr::pull(symbol) %>% unique()
  if (length(genes_all) < 5) { warning("Too few genes after |LFC|: ", coef_name); return(invisible(NULL)) }
  
  ego <- enrichGO(gene = genes_all, universe = universe_syms, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                  ont = GO$ont_choice, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                  minGSSize = GO$minGSSize, maxGSSize = GO$maxGSSize, readable = TRUE)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) { warning("No GO terms: ", coef_name); return(invisible(NULL)) }
  
  ego_s <- if (GO$do_simplify) {
    tryCatch(simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang"),
             error = function(e) { warning("simplify() failed: ", e$message); ego })
  } else ego
  
  base <- safe_name(coef_name)
  write.csv(as.data.frame(ego),   file.path(outs$go_tables_raw,   paste0(base, ".csv")), row.names = FALSE)
  write.csv(as.data.frame(ego_s), file.path(outs$go_tables_simpl, paste0(base, ".csv")), row.names = FALSE)
  
  p_bar_raw <- barplot(ego, showCategory = 10, drop = TRUE, font.size = 4, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y")
  p_bar_raw <- wrap_go_labels(p_bar_raw)
  p_bar_raw <- add_title(p_bar_raw, paste0("GO (", GO$ont_choice, ") — bar (raw): ", coef_name))
  save_plot_both(file.path(outs$go_plots_bar_raw, paste0(base, "_bar_raw")), p_bar_raw, w = GO$W, h = GO$H, dpi = GO$DPI)
  
  p_dot_raw <- dotplot(ego, showCategory = 10, split = "ONTOLOGY", font.size = 4) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y")
  p_dot_raw <- wrap_go_labels(p_dot_raw)
  p_dot_raw <- shrink_dots(p_dot_raw, range = c(1.2, 4))
  p_dot_raw <- add_title(p_dot_raw, paste0("GO (", GO$ont_choice, ") — dot (raw): ", coef_name))
  save_plot_both(file.path(outs$go_plots_dot_raw, paste0(base, "_dot_raw")), p_dot_raw, w = GO$W, h = GO$H, dpi = GO$DPI)
  
  p_bar_s <- barplot(ego_s, showCategory = 10, drop = TRUE, font.size = 4, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y")
  p_bar_s <- wrap_go_labels(p_bar_s)
  p_bar_s <- add_title(p_bar_s, paste0("GO (", GO$ont_choice, ") — bar (simplified): ", coef_name))
  save_plot_both(file.path(outs$go_plots_bar_simpl, paste0(base, "_bar_simpl")), p_bar_s, w = GO$W, h = GO$H, dpi = GO$DPI)
  
  p_dot_s <- dotplot(ego_s, showCategory = 10, split = "ONTOLOGY", font.size = 4) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y")
  p_dot_s <- wrap_go_labels(p_dot_s)
  p_dot_s <- shrink_dots(p_dot_s, range = c(1.2, 4))
  p_dot_s <- add_title(p_dot_s, paste0("GO (", GO$ont_choice, ") — dot (simplified): ", coef_name))
  save_plot_both(file.path(outs$go_plots_dot_simpl, paste0(base, "_dot_simpl")), p_dot_s, w = GO$W, h = GO$H, dpi = GO$DPI)
  
  ab <- parse_AB(coef_name, main_var); A <- ab$A; B <- ab$B
  sig_up_A <- sym_map %>% dplyr::filter(lfc >  thr$go_lfc_cut) %>% dplyr::pull(symbol)
  sig_up_B <- sym_map %>% dplyr::filter(lfc < -thr$go_lfc_cut) %>% dplyr::pull(symbol)
  
  gene_sets <- list()
  if (length(sig_up_A) >= GO$min_genes_dir) gene_sets[[paste0("Up_in_", A)]] <- sig_up_A
  if (length(sig_up_B) >= GO$min_genes_dir) gene_sets[[paste0("Up_in_", B)]] <- sig_up_B
  
  if (length(gene_sets)) {
    cc <- compareCluster(geneCluster = gene_sets, fun = "enrichGO", OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                         ont = GO$ont_choice, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                         minGSSize = GO$minGSSize, maxGSSize = GO$maxGSSize, readable = TRUE)
    write.csv(as.data.frame(cc), file.path(outs$go_tables_cc, paste0(base, "_compareCluster.csv")), row.names = FALSE)
    
    p_cc <- dotplot(cc, showCategory = 10, split = "ONTOLOGY", font.size = 3) +
      facet_grid(ONTOLOGY ~ ., scales = "free_y")
    p_cc <- wrap_go_labels(p_cc)
    p_cc <- shrink_dots(p_cc, range = c(1.0, 3.5))
    p_cc <- add_title(p_cc, paste0("GO (", GO$ont_choice, ") — ", A, " vs ", B),
                      "Direction of enrichment: Up_in_<group>")
    save_plot_both(file.path(outs$go_plots_cc, paste0(base, "_compareCluster_dot")), p_cc, w = GO$W, h = GO$H, dpi = GO$DPI)
  }
  
  invisible(TRUE)
}

##### FIND COEFFICIENT NAMES FROM RESULTS #####
coef_names <- grep(paste0("^", main_var, "_.*_vs_"), resultsNames(dds_grouped), value = TRUE)
message("Found ", length(coef_names), " contrasts for ", main_var, ".")
invisible(lapply(coef_names, run_go_for_coef))
message("\nAll done. Outputs in: ", normalizePath(paths$out_base, mustWork = FALSE))