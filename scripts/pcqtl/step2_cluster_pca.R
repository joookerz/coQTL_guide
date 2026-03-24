#!/usr/bin/env Rscript
# ============================================================================
# Step 2: Single-cell PCA per gene cluster
# ============================================================================
# Optimisation: only loads genes that appear in clusters + covariate columns
# from the readcount file, dramatically reducing memory for large cell types.
# Output per cluster: pheno_with_pcs.tsv  (individual + covariates + PC scores)
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

cat(sprintf("\n=== Cluster PCA (Single Cell) ===\nStart: %s\n\n", Sys.time()))

# ── config ──────────────────────────────────────────────────────────────────
get_this_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("--file=", cmd_args)
  if (length(file_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[file_idx])))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  normalizePath(".")
}

script_dir <- dirname(get_this_path())
source(file.path(script_dir, "load_config.R"))
cfg <- load_pcqtl_config(script_dir)

STEP2_DIR          <- file.path(cfg$pcqtl_dir, "step2_pca")
dir.create(STEP2_DIR, showWarnings = FALSE, recursive = TRUE)

COVAR_COLS         <- cfg$covar_cols
VARIANCE_THRESHOLD <- cfg$variance_threshold

# ── cluster info ────────────────────────────────────────────────────────────
cat("Loading cluster info ...\n")
clusters      <- fread(cfg$cluster_file)
cluster_genes <- fread(cfg$cluster_gene_file)
cat(sprintf("  Clusters: %d | Unique genes in clusters: %d\n\n",
            nrow(clusters), length(unique(cluster_genes$gene_name))))

# ── column selection (memory-optimised read) ───────────────────────────────
cat("Reading column names from readcount file ...\n")
header_cols   <- names(fread(cfg$count_file, nrows = 0))
all_genes     <- unique(cluster_genes$gene_name)
genes_avail   <- intersect(all_genes,  header_cols)
covars_avail  <- intersect(COVAR_COLS, header_cols)
individual_col <- intersect(cfg$individual_col_candidates, header_cols)
if (!length(individual_col)) {
  stop("Count file must contain one individual column: ", paste(cfg$individual_col_candidates, collapse = ", "))
}
individual_col <- individual_col[1]
load_cols     <- unique(c(individual_col, genes_avail, covars_avail))
cat(sprintf("  Selecting %d genes + %d covariate cols (total file cols: %d)\n\n",
            length(genes_avail), length(covars_avail), length(header_cols)))

# ── load data ───────────────────────────────────────────────────────────────
cat("Loading readcount data (selected columns only) ...\n")
t0        <- Sys.time()
expr_data <- fread(cfg$count_file, select = load_cols)
cat(sprintf("  Loaded in %.2f min — %d cells × %d columns\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            nrow(expr_data), ncol(expr_data)))

sample_ids <- expr_data[[individual_col]]

# ── per-cluster PCA ─────────────────────────────────────────────────────────
cat("Running PCA per cluster ...\n")

process_cluster <- function(i, crow) {
  cid        <- crow$cluster_id
  cgenes     <- strsplit(crow$genes, ",")[[1]]
  cluster_dir <- file.path(STEP2_DIR, cid)
  dir.create(cluster_dir, showWarnings = FALSE, recursive = TRUE)

  avail <- cgenes[cgenes %in% colnames(expr_data)]
  if (length(avail) < 2) {
    return(data.table(cluster_id = cid, status = "FAILED",
                      reason = "Too few genes available", n_genes = length(avail)))
  }

  cluster_expr <- as.matrix(expr_data[, ..avail])
  keep         <- complete.cases(cluster_expr)
  cluster_expr <- cluster_expr[keep, , drop = FALSE]

  n_samples <- nrow(cluster_expr)
  n_genes   <- ncol(cluster_expr)
  if (n_samples < n_genes) {
    return(data.table(cluster_id = cid, status = "FAILED",
                      reason = "n_samples < n_genes", n_genes = n_genes, n_samples = n_samples))
  }

  tryCatch({
    pca <- prcomp(cluster_expr, center = TRUE, scale. = FALSE)

    var_expl    <- pca$sdev^2 / sum(pca$sdev^2)
    cum_var     <- cumsum(var_expl)
    n_pcs_95    <- which(cum_var >= VARIANCE_THRESHOLD)[1L]
    if (is.na(n_pcs_95)) n_pcs_95 <- length(var_expl)
    n_pcs_total <- ncol(pca$x)

    # save PCA object
    save(pca, file = file.path(cluster_dir, "pca_results.rda"))

    # phenotype file: individual + covariates + all PCs
    pheno <- data.table(
      individual = sample_ids[keep],
      expr_data[keep, ..covars_avail],
      as.data.table(pca$x)          # columns PC1, PC2, …  (uppercase)
    )
    fwrite(pheno, file.path(cluster_dir, "pheno_with_pcs.tsv"), sep = "\t", quote = FALSE)

    # variance explained
    fwrite(
      data.table(PC           = paste0("PC", seq_len(n_pcs_total)),
                 Eigenvalue   = pca$sdev^2,
                 VarExplained = var_expl,
                 CumVar       = cum_var),
      file.path(cluster_dir, "variance_explained.tsv"), sep = "\t", quote = FALSE)

    # gene loadings
    fwrite(
      data.table(Gene = rownames(pca$rotation), as.data.table(pca$rotation)),
      file.path(cluster_dir, "gene_loadings.tsv"), sep = "\t", quote = FALSE)

    data.table(cluster_id   = cid,
               status       = "SUCCESS",
               n_genes      = n_genes,
               n_samples    = n_samples,
               n_pcs_total  = n_pcs_total,
               n_pcs_95pct  = n_pcs_95,
               var_95pct    = cum_var[n_pcs_95])
  }, error = function(e) {
    data.table(cluster_id = cid, status = "FAILED", reason = as.character(e),
               n_genes = n_genes, n_samples = n_samples)
  })
}

results <- rbindlist(lapply(seq_len(nrow(clusters)), function(i) {
  if (i %% 100 == 0 || i == nrow(clusters))
    cat(sprintf("  %d / %d clusters (%.0f%%)\n", i, nrow(clusters), 100 * i / nrow(clusters)))
  process_cluster(i, clusters[i])
}))

# ── summary ─────────────────────────────────────────────────────────────────
n_ok   <- sum(results$status == "SUCCESS")
n_fail <- sum(results$status == "FAILED")
cat(sprintf("\nResults: %d SUCCESS, %d FAILED\n", n_ok, n_fail))
if (n_fail > 0) {
  cat("Failed:\n")
  print(results[status == "FAILED", .(cluster_id, reason)])
}
fwrite(results, file.path(STEP2_DIR, "all_clusters_summary.tsv"), sep = "\t", quote = FALSE)

# ── PC list for SAIGE ───────────────────────────────────────────────────────
pc_list <- rbindlist(lapply(which(results$status == "SUCCESS"), function(i) {
  cid  <- results$cluster_id[i]
  npcs <- results$n_pcs_95pct[i]
  data.table(cluster_id = cid,
             PC         = paste0("PC", seq_len(npcs)),
             pheno_file = file.path(STEP2_DIR, cid, "pheno_with_pcs.tsv"))
}))
fwrite(pc_list, file.path(STEP2_DIR, "pc_list_for_saige.tsv"), sep = "\t", quote = FALSE)

cat(sprintf("\nTotal cluster-PCs for SAIGE: %d\n", nrow(pc_list)))
cat(sprintf("Output dir : %s\n", STEP2_DIR))
cat(sprintf("End        : %s\n", Sys.time()))
