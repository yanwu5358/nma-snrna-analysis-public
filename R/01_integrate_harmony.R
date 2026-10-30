#!/usr/bin/env Rscript
# =============================================================================
# 01_integrate_harmony.R
#
# Stage 2 of the analysis pipeline: merge the 9 cleaned per-sample objects,
# normalize with SCTransform v2 (parallel per sample), integrate with Harmony
# on batch, and produce a clustered UMAP.
#
# Example placeholders only: replace with your own sample metadata and file names.
# Design: 3 batches × 3 conditions
#   Sample ID         batch   condition
#   sampleB1_101         1      CTRL
#   sampleB1_102         1      FUS
#   sampleB1_103         1      TDP
#   sampleB2_101         2      CTRL
#   ...                ...    ...
#   sampleB3_103         3      TDP
#
#
# Usage:
#   Rscript 01_integrate_harmony.R \
#     --rds-root ./data/cleaned \
#     --out-dir  ./results/integration
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(harmony)
  library(future)
  library(future.apply)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

option_list <- list(
  make_option("--rds-root", type = "character",
              default = "./data/cleaned",
              help = "Root dir with {sample}/rds/{sample}_seu_clean.rds files"),
  make_option("--out-dir",  type = "character",
              default = "./results/integration",
              help = "Output root directory (a subfolder per group choice is created)"),
  make_option("--group-by-vars", type = "character", default = "batch",
              help = paste("Comma-separated metadata column(s) for Harmony to",
                           "correct on. 'batch' (gentle), 'sample' (aggressive),",
                           "'batch,sample' (nested). [default: %default]")),
  make_option("--theta",    type = "double",    default = 2.0,
              help = paste("Harmony theta — correction strength.",
                           "Higher = stronger mixing. Lower for sample-level",
                           "correction when sample is confounded with condition.",
                           "[default: %default]")),
  make_option("--workers",  type = "integer",   default = 9L,
              help = "Parallel SCT workers [default: %default]"),
  make_option("--nfeatures",type = "integer",   default = 3000L,
              help = "SelectIntegrationFeatures nfeatures [default: %default]"),
  make_option("--ndims",    type = "integer",   default = 30L,
              help = "PCs used for Harmony/UMAP [default: %default]"),
  make_option("--res",      type = "double",    default = 0.5,
              help = "Cluster resolution [default: %default]"),
  make_option("--force",    action = "store_true", default = FALSE,
              help = "Recompute SCT+merge+PCA even if checkpoint exists")
)
opts <- parse_args(OptionParser(option_list = option_list))
opts$rds_root <- normalizePath(path.expand(opts$`rds-root`),  mustWork = FALSE)
opts$out_dir  <- normalizePath(path.expand(opts$`out-dir`),   mustWork = FALSE)

# Parse comma-separated group-by-vars into a character vector
group_vars <- trimws(strsplit(opts$`group-by-vars`, ",", fixed = TRUE)[[1]])
allowed_vars <- c("batch", "sample", "condition")
bad <- setdiff(group_vars, allowed_vars)
if (length(bad)) stop("Unknown group var(s): ", paste(bad, collapse = ", "),
                      "  (allowed: batch, sample, condition)")
if ("condition" %in% group_vars) {
  warning("Correcting on 'condition' will erase the very signal you want to compare. ",
          "Are you sure?")
}
# Tag for subfolder naming, e.g. "by_batch" or "by_batch_sample"
group_tag <- paste0("by_", paste(group_vars, collapse = "_"))

# Each run writes to its own subfolder under out_dir
run_dir  <- file.path(opts$out_dir, group_tag)
plot_dir <- file.path(run_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
rds_out  <- file.path(run_dir, "integrated_harmony.rds")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
              paste0(..., collapse = "")))
  flush.console()
}

t0 <- Sys.time()

# -----------------------------------------------------------------------------
# 1. Sample manifest (parsed from sample IDs)
# -----------------------------------------------------------------------------
SAMPLES <- c("sampleB1_101", "sampleB1_102", "sampleB1_103",
             "sampleB2_101", "sampleB2_102", "sampleB2_103",
             "sampleB3_101", "sampleB3_102", "sampleB3_103")

cond_map <- c("101" = "CTRL", "102" = "FUS", "103" = "TDP")

manifest <- data.frame(
  sample = SAMPLES,
  batch  = sub("^sampleB([0-9]+)_.*$", "\\1", SAMPLES),
  rep    = sub("^.*_",            "",     SAMPLES),
  stringsAsFactors = FALSE
)
manifest$condition <- factor(cond_map[manifest$rep],
                             levels = c("CTRL", "FUS", "TDP"))
manifest$rds_path <- file.path(opts$rds_root, manifest$sample, "rds",
                               paste0(manifest$sample, "_seu_clean.rds"))

log_msg("Sample manifest:")
print(manifest[, c("sample", "batch", "condition", "rep")])

# Sanity check: every input file present
missing_files <- manifest$rds_path[!file.exists(manifest$rds_path)]
if (length(missing_files)) {
  stop("Missing cleaned RDS files:\n  ", paste(missing_files, collapse = "\n  "))
}

# -----------------------------------------------------------------------------
# 2. Parallel SCTransform per sample
# -----------------------------------------------------------------------------
# Each worker is a separate R session — isolated memory, no copy-on-modify issues.
# future.globals.maxSize bumped to 50GB to handle large per-sample matrices.
log_msg(sprintf("Setting up %d-worker future plan ...", opts$workers))
plan(multisession, workers = opts$workers)
options(future.globals.maxSize = 50 * 1024^3)

sct_one <- function(i, manifest) {
  suppressPackageStartupMessages({
    library(Seurat); library(Matrix)
  })
  row <- manifest[i, ]
  seu <- readRDS(row$rds_path)

  # Tag metadata
  seu$sample    <- row$sample
  seu$batch     <- row$batch
  seu$condition <- row$condition
  seu$rep       <- row$rep

  # Drop downstream-derived columns we'll recompute on the integrated object
  for (col in c("seurat_clusters", "SCT_snn_res.0.5", "CellBender_snn_res.0.5")) {
    if (col %in% colnames(seu@meta.data)) seu[[col]] <- NULL
  }

  # CellBender counts can be fractional; SCT wants non-negative integers.
  # Round in place — typical loss is <1% per cell, well below technical noise.
  DefaultAssay(seu) <- "CellBender"
  cts <- GetAssayData(seu, layer = "counts", assay = "CellBender")
  cts@x <- round(cts@x)
  seu <- SetAssayData(seu, layer = "counts", assay = "CellBender", new.data = cts)

  # SCTransform v2 — depth-aware normalization. Suppress the (chatty) progress.
  seu <- SCTransform(seu, assay = "CellBender", vst.flavor = "v2",
                     return.only.var.genes = FALSE, verbose = FALSE)

  # Keep only what we need downstream to shrink the in-memory footprint
  DefaultAssay(seu) <- "SCT"
  DietSeurat(seu, assays = c("SCT", "CellBender"), dimreducs = NULL, graphs = NULL)
}

# -----------------------------------------------------------------------------
# 3. Merge + select integration features + PCA   (with a shared CHECKPOINT)
# -----------------------------------------------------------------------------
# SCT + merge + PCA are the expensive, deterministic, group_var-independent
# stages. We cache the post-PCA object at the parent out_dir level so reruns
# with different Harmony params skip straight to step 4.
ckpt_rds <- file.path(opts$out_dir, "merged_sct_pca.rds")

if (file.exists(ckpt_rds) && !isTRUE(opts$force)) {
  log_msg("Loading cached post-PCA object: ", ckpt_rds,
          "  (delete this file or pass --force to recompute SCT+PCA)")
  seu_m <- readRDS(ckpt_rds)
} else {
  log_msg(sprintf("Running SCTransform per sample (%d in parallel) ...",
                  opts$workers))
  seu_list <- future_lapply(seq_len(nrow(manifest)),
                            sct_one, manifest = manifest,
                            future.seed = TRUE)
  names(seu_list) <- manifest$sample
  log_msg("All per-sample SCT done.")

  # Drop the parallel plan before doing heavy single-process work
  plan(sequential); gc(verbose = FALSE)

  log_msg("Selecting integration features ...")
  features <- SelectIntegrationFeatures(object.list = seu_list,
                                        nfeatures = opts$nfeatures)
  log_msg(sprintf("  %d features selected.", length(features)))

  log_msg("Merging objects ...")
  seu_m <- merge(seu_list[[1]], seu_list[-1],
                 add.cell.ids = names(seu_list),
                 project = "NMA_integrated")
  rm(seu_list); gc(verbose = FALSE)

  DefaultAssay(seu_m) <- "SCT"
  VariableFeatures(seu_m) <- features

  log_msg("Running PCA ...")
  seu_m <- RunPCA(seu_m, features = features, npcs = 50, verbose = FALSE)

  log_msg("Saving merged+SCT+PCA checkpoint to ", ckpt_rds)
  saveRDS(seu_m, ckpt_rds)
}

# Always (re-)write the elbow plot into the per-run folder
ggsave(file.path(plot_dir, "elbow_PCA.pdf"),
       ElbowPlot(seu_m, ndims = 50) +
         geom_vline(xintercept = opts$ndims, linetype = "dashed", colour = "red"),
       width = 7, height = 5)

# -----------------------------------------------------------------------------
# 4. Harmony integration
# -----------------------------------------------------------------------------
log_msg(sprintf("Running Harmony on %s (theta=%g, %d PCs) ...",
                paste(group_vars, collapse = "+"), opts$theta, opts$ndims))
seu_m <- RunHarmony(seu_m,
                    group.by.vars = group_vars,
                    theta         = opts$theta,
                    dims.use      = 1:opts$ndims,
                    plot_convergence = FALSE,
                    verbose       = FALSE)

# -----------------------------------------------------------------------------
# 5. UMAP + clustering on the Harmony embedding
# -----------------------------------------------------------------------------
log_msg("UMAP + clustering on Harmony reduction ...")
seu_m <- RunUMAP(seu_m,        reduction = "harmony", dims = 1:opts$ndims,
                 verbose = FALSE)
seu_m <- FindNeighbors(seu_m,  reduction = "harmony", dims = 1:opts$ndims,
                       verbose = FALSE)
seu_m <- FindClusters(seu_m,   resolution = opts$res, verbose = FALSE)

# -----------------------------------------------------------------------------
# 6. QC plots — full before/after across BOTCH technical and biological vars
# -----------------------------------------------------------------------------
log_msg("Writing integration QC plots ...")
DefaultAssay(seu_m) <- "SCT"

# Always show PCA-vs-Harmony mixing for batch, sample, AND condition so the
# two runs (by_batch / by_sample) are directly comparable.
diag_vars <- list(
  batch     = NULL,
  sample    = NULL,
  condition = c(CTRL = "grey50", FUS = "#E41A1C", TDP = "#377EB8")
)
mixing_panels <- list()
for (gv in names(diag_vars)) {
  p_pca <- DimPlot(seu_m, reduction = "pca",     group.by = gv, pt.size = 0.1) +
           ggtitle(paste0("PCA · by ", gv))
  p_har <- DimPlot(seu_m, reduction = "harmony", group.by = gv, pt.size = 0.1) +
           ggtitle(paste0("Harmony · by ", gv))
  if (!is.null(diag_vars[[gv]])) {
    p_pca <- p_pca + scale_colour_manual(values = diag_vars[[gv]])
    p_har <- p_har + scale_colour_manual(values = diag_vars[[gv]])
  }
  mixing_panels[[gv]] <- p_pca | p_har
}
ggsave(file.path(plot_dir, "harmony_mixing_full.pdf"),
       wrap_plots(mixing_panels, ncol = 1) +
         plot_annotation(title = paste0("Harmony group.by.vars = ",
                                        paste(group_vars, collapse = "+"),
                                        ",  theta = ", opts$theta)),
       width = 14, height = 18)

# UMAP coloured by each metadata of interest
make_umap <- function(grp, title = grp, cols = NULL) {
  p <- DimPlot(seu_m, reduction = "umap", group.by = grp, pt.size = 0.1) +
       ggtitle(title)
  if (!is.null(cols)) p <- p + scale_colour_manual(values = cols)
  p
}
cond_cols <- c(CTRL = "grey50", FUS = "#E41A1C", TDP = "#377EB8")

ggsave(file.path(plot_dir, "umap_overview.pdf"),
       (make_umap("seurat_clusters", "Clusters") |
        make_umap("batch",          "Batch")) /
       (make_umap("condition",      "Condition", cond_cols) |
        make_umap("sample",         "Sample")),
       width = 16, height = 12)

# Per-condition UMAP (split view)
ggsave(file.path(plot_dir, "umap_split_by_condition.pdf"),
       DimPlot(seu_m, reduction = "umap", group.by = "seurat_clusters",
               split.by = "condition", pt.size = 0.1, ncol = 3) +
         ggtitle("Clusters · split by condition"),
       width = 18, height = 6)

# Cluster composition by condition (proportions)
comp <- seu_m@meta.data |>
  dplyr::count(seurat_clusters, condition) |>
  dplyr::group_by(seurat_clusters) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::ungroup()
ggsave(file.path(plot_dir, "cluster_composition_by_condition.pdf"),
       ggplot(comp, aes(seurat_clusters, prop, fill = condition)) +
         geom_col() +
         scale_fill_manual(values = cond_cols) +
         labs(x = "Cluster", y = "Proportion",
              title = "Cluster composition by condition") +
         theme_classic(),
       width = 10, height = 5)

# Write a summary table for quick review
write.csv(comp, file.path(run_dir, "cluster_composition.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Save the integrated object
# -----------------------------------------------------------------------------
log_msg("Saving integrated object ...")
saveRDS(seu_m, rds_out)

dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
log_msg(sprintf("DONE in %.1f min.  Object saved to %s", dt, rds_out))
log_msg(sprintf("  %d cells across %d clusters",
                ncol(seu_m), length(unique(seu_m$seurat_clusters))))
