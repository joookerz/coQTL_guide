#!/usr/bin/env Rscript

# =============================================================================
# Method 2 - Step 2: Calculate SC hurdle associations (CHUNKED VERSION)
# =============================================================================
# Process a subset of genes on one chromosome to keep runtime under 6 hours
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(fasthurdle)
})

require_fasthurdle_version <- function(min_version = "1.1.1") {
  installed <- tryCatch(utils::packageVersion("fasthurdle"), error = function(e) NULL)
  if (is.null(installed)) {
    stop("Package 'fasthurdle' is not installed.")
  }
  if (installed < package_version(min_version)) {
    stop(sprintf(
      "fasthurdle >= %s is required (found %s). Please update fasthurdle.",
      min_version, as.character(installed)
    ))
  }
  as.character(installed)
}

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
cfg <- load_config(script_dir)
fasthurdle_version <- require_fasthurdle_version("1.1.1")

option_list <- list(
  make_option("--chr", type = "integer", help = "Chromosome number"),
  make_option("--gene_start", type = "integer", help = "Start gene index (1-based)"),
  make_option("--gene_end", type = "integer", help = "End gene index (inclusive)"),
  make_option("--chunk_id", type = "integer", default = 1, help = "Chunk ID for naming"),
  make_option("--p_threshold", type = "double", default = NA_real_),
  make_option("--nonzero_cutoff", type = "double", default = NA_real_)
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$chr) || is.null(opt$gene_start) || is.null(opt$gene_end)) {
  stop("Please provide --chr, --gene_start, and --gene_end")
}

CHR <- opt$chr
GENE_START <- opt$gene_start
GENE_END <- opt$gene_end
CHUNK_ID <- opt$chunk_id
P_THRESHOLD <- if (is.na(opt$p_threshold)) cfg$p_threshold else opt$p_threshold
NONZERO_CUTOFF <- if (is.na(opt$nonzero_cutoff)) cfg$nonzero_cutoff else opt$nonzero_cutoff

cat(sprintf("\n=== Method 2 Step 2: SC Hurdle (Chunked) - Chr%d Chunk%d (%s) ===\n",
            CHR, CHUNK_ID, cfg$cell_type))
cat(sprintf("fasthurdle version: %s\n", fasthurdle_version))
cat(sprintf("Gene range: %d-%d\n", GENE_START, GENE_END))
cat(sprintf("Start time: %s\n\n", Sys.time()))

COUNT_FILE <- cfg$count_file
GENE_INFO_FILE <- cfg$gene_info_file
FILTER_FILE <- file.path(cfg$method2_results, "filtered_genes.tsv")

ASSOC_DIR <- cfg$chunk_dir
dir.create(ASSOC_DIR, recursive = TRUE, showWarnings = FALSE)

# Load gene metadata
cat("Loading gene metadata...\n")
gene_info <- fread(GENE_INFO_FILE)
filt <- fread(FILTER_FILE)
merged <- merge(gene_info, filt[keep == TRUE, .(gene_name, nonzero_prop)], by = "gene_name")
merged <- merged[!duplicated(gene_name)]
global_counts <- merged[, .N, by = chr_numeric]
global_counts <- global_counts[chr_numeric %in% seq_len(22)]
global_pairs <- global_counts[, sum(N * (N - 1) / 2)]
bonferroni_threshold <- P_THRESHOLD / global_pairs

chr_genes_all <- merged[chr_numeric == CHR][order(start), gene_name]

n_genes_chr <- length(chr_genes_all)
cat(sprintf("  Chromosome %d total filtered genes: %d\n", CHR, n_genes_chr))

# Select chunk
if (GENE_END > n_genes_chr) GENE_END <- n_genes_chr
chr_genes <- chr_genes_all[GENE_START:GENE_END]
n_genes <- length(chr_genes)

cat(sprintf("  Selected genes %d-%d: %d genes\n", GENE_START, GENE_END, n_genes))

if (n_genes < 2) {
  cat("  Not enough genes; skipping.\n")
  quit(status = 0)
}

# Load counts
cat("Loading single-cell counts...\n")
header_cols <- names(fread(COUNT_FILE, nrows = 0))
id_col <- intersect(c("CellID", "barcode"), header_cols)
if (!length(id_col)) stop("Count file must contain CellID or barcode column")
select_cols <- unique(c(id_col[1], chr_genes))
count_dt <- fread(COUNT_FILE, select = select_cols)
count_mat <- as.matrix(count_dt[, ..chr_genes])
n_cells <- nrow(count_mat)

# Calculate associations
pair_idx <- combn(n_genes, 2, simplify = FALSE)
n_pairs <- length(pair_idx)
cat(sprintf("  Cells: %d, Genes: %d, Pairs: %d\n", n_cells, n_genes, n_pairs))

run_pair <- function(i, j) {
  df <- data.frame(count_i = count_mat[, i], count_j = count_mat[, j])
  if (all(df$count_i == 0) || all(df$count_j == 0)) {
    return(list(p_count = NA_real_, p_zero = NA_real_))
  }

  fit_ij <- tryCatch(fasthurdle(count_i ~ count_j, data = df,
                                dist = "poisson", zero.dist = "binomial"),
                     error = function(e) NULL)
  fit_ji <- tryCatch(fasthurdle(count_j ~ count_i, data = df,
                                dist = "poisson", zero.dist = "binomial"),
                     error = function(e) NULL)

  if (is.null(fit_ij) && is.null(fit_ji)) {
    return(list(p_count = NA_real_, p_zero = NA_real_))
  }

  grab <- function(model, component, term) {
    if (is.null(model)) return(NA_real_)
    smry <- summary(model)$coefficients
    if (!is.null(smry[[component]]) && term %in% rownames(smry[[component]])) {
      return(smry[[component]][term, "Pr(>|z|)"])
    }
    NA_real_
  }

  p_count_vals <- c(grab(fit_ij, "count", "count_j"), grab(fit_ji, "count", "count_i"))
  p_zero_vals <- c(grab(fit_ij, "zero", "count_j"), grab(fit_ji, "zero", "count_i"))

  p_count <- min(p_count_vals, na.rm = TRUE)
  if (!is.finite(p_count)) p_count <- NA_real_
  p_zero <- min(p_zero_vals, na.rm = TRUE)
  if (!is.finite(p_zero)) p_zero <- NA_real_

  list(p_count = p_count, p_zero = p_zero)
}

pval_count <- matrix(NA_real_, n_genes, n_genes, dimnames = list(chr_genes, chr_genes))
pval_zero <- matrix(NA_real_, n_genes, n_genes, dimnames = list(chr_genes, chr_genes))

cat("Calculating associations...\n")
start_time <- Sys.time()

for (idx in seq_along(pair_idx)) {
  pair <- pair_idx[[idx]]
  i <- pair[1]; j <- pair[2]
  res <- run_pair(i, j)
  pval_count[i, j] <- pval_count[j, i] <- res$p_count
  pval_zero[i, j] <- pval_zero[j, i] <- res$p_zero

  if (idx %% 1000 == 0 || idx == n_pairs) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- idx / max(elapsed, 1e-6)
    eta <- (n_pairs - idx) / max(rate, 1e-6)
    cat(sprintf("  %d/%d pairs (%.1f%%) -- %.2f pairs/s, ETA %.1f min\n",
                idx, n_pairs, 100 * idx / n_pairs, rate, eta / 60))
  }
}

# Extract significant pairs
sig_pairs <- data.table()
for (i in 1:(n_genes - 1)) {
  for (j in (i + 1):n_genes) {
    p_count <- pval_count[i, j]
    p_zero <- pval_zero[i, j]
    # Check for NA before comparison to avoid "missing value where TRUE/FALSE needed"
    if ((!is.na(p_count) && p_count < bonferroni_threshold) ||
        (!is.na(p_zero) && p_zero < bonferroni_threshold)) {
      sig_pairs <- rbind(sig_pairs, data.table(
        Gene1 = chr_genes[i],
        Gene2 = chr_genes[j],
        Pvalue_count = p_count,
        Pvalue_zero = p_zero
      ))
    }
  }
}

# Save results
out_base <- file.path(ASSOC_DIR, sprintf("chr%d_chunk%03d", CHR, CHUNK_ID))
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

fwrite(as.data.table(pval_count, keep.rownames = "Gene"),
       file.path(out_base, "pvalues_count.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(as.data.table(pval_zero, keep.rownames = "Gene"),
       file.path(out_base, "pvalues_zero.tsv.gz"), sep = "\t", compress = "gzip")

if (nrow(sig_pairs) > 0) {
  fwrite(sig_pairs, file.path(out_base, "significant_pairs.tsv"), sep = "\t")
}

summary_dt <- data.table(
  Chromosome = CHR,
  ChunkID = CHUNK_ID,
  GeneStart = GENE_START,
  GeneEnd = GENE_END,
  NumGenes = n_genes,
  NumPairs = n_pairs,
  GlobalPairs = global_pairs,
  BonferroniThreshold = bonferroni_threshold,
  SignificantPairs = nrow(sig_pairs)
)
fwrite(summary_dt, file.path(out_base, "summary.tsv"), sep = "\t")

cat("\n=== Complete ===\n")
cat(sprintf("End time: %s\n", Sys.time()))
cat(sprintf("Output: %s\n", out_base))
