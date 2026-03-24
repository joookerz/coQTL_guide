#!/usr/bin/env Rscript

# Reapply global Bonferroni threshold using existing p-value matrices (no refitting)

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

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

option_list <- list(
  make_option("--chunk_dir", default = NULL,
              help = "Directory containing chunk outputs with pvalues_count/pvalues_zero"),
  make_option("--filter_file", default = NULL,
              help = "Filtered genes file from Step 1"),
  make_option("--gene_info", default = NULL,
              help = "Gene information table"),
  make_option("--alpha", type = "double", default = NA_real_,
              help = "Family-wise alpha for Bonferroni")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$chunk_dir)) opt$chunk_dir <- cfg$chunk_dir
if (is.null(opt$filter_file)) opt$filter_file <- file.path(cfg$method2_results, "filtered_genes.tsv")
if (is.null(opt$gene_info)) opt$gene_info <- cfg$gene_info_file
if (is.na(opt$alpha)) opt$alpha <- cfg$p_threshold

cat(sprintf("\n=== Reapplying SC hurdle Bonferroni threshold (%s) ===\n",
            cfg$cell_type))
cat(sprintf("Chunk dir: %s\n", opt$chunk_dir))

if (!dir.exists(opt$chunk_dir)) {
  stop("Chunk directory not found. Run chunked SC hurdle first.")
}

if (!file.exists(opt$filter_file)) {
  stop("Filtered genes file not found. Run step1_filter_sparse_genes.R first.")
}

gene_info <- fread(opt$gene_info)
filt <- fread(opt$filter_file)[keep == TRUE, .(gene_name)]
merged <- merge(gene_info, filt, by = "gene_name")
merged <- merged[!duplicated(gene_name)]
global_counts <- merged[, .N, by = chr_numeric]
global_counts <- global_counts[chr_numeric %in% seq_len(22)]
global_pairs <- global_counts[, sum(N * (N - 1) / 2)]
bonf_thr <- opt$alpha / global_pairs

cat(sprintf("Total dense genes: %d  => Global pairs: %s\n",
            nrow(merged),
            format(global_pairs, big.mark = ",")))
cat(sprintf("Bonferroni threshold: %.6g\n\n", bonf_thr))

chunk_dirs <- sort(list.dirs(opt$chunk_dir, recursive = FALSE, full.names = TRUE))
if (!length(chunk_dirs)) {
  stop("No chunk subdirectories detected.")
}

process_chunk <- function(chunk_path) {
  count_file <- file.path(chunk_path, "pvalues_count.tsv.gz")
  zero_file <- file.path(chunk_path, "pvalues_zero.tsv.gz")
  if (!file.exists(count_file) || !file.exists(zero_file)) {
    cat("  Skipping chunk (missing pvalue files):", basename(chunk_path), "\n")
    return(NULL)
  }
  count_dt <- fread(cmd = sprintf("gzip -dc %s", count_file))
  zero_dt <- fread(cmd = sprintf("gzip -dc %s", zero_file))
  genes <- count_dt[[1]]
  count_mat <- as.matrix(count_dt[, -1])
  zero_mat <- as.matrix(zero_dt[, -1])
  rownames(count_mat) <- genes
  colnames(count_mat) <- genes
  rownames(zero_mat) <- genes
  colnames(zero_mat) <- genes

  sig_idx <- which((count_mat < bonf_thr) | (zero_mat < bonf_thr), arr.ind = TRUE)
  sig_idx <- sig_idx[sig_idx[,1] < sig_idx[,2], , drop = FALSE]

  if (nrow(sig_idx) == 0) {
    fwrite(data.table(), file.path(chunk_path, "significant_pairs.tsv"), sep = "\t")
    fwrite(data.table(
             Chromosome = sub(".*chr(\\d+)_chunk.*", "\\1", basename(chunk_path)),
             ChunkID = sub(".*chunk(\\d+)", "\\1", basename(chunk_path)),
             SignificantPairs = 0,
             BonferroniThreshold = bonf_thr),
           file.path(chunk_path, "summary.tsv"), sep = "\t")
    return(0L)
  }

  pairs_dt <- data.table(
    Gene1 = genes[sig_idx[,1]],
    Gene2 = genes[sig_idx[,2]],
    Pvalue_count = count_mat[sig_idx],
    Pvalue_zero = zero_mat[sig_idx]
  )
  fwrite(pairs_dt, file.path(chunk_path, "significant_pairs.tsv"), sep = "\t")

  summ_dt <- data.table(
    Chromosome = sub(".*chr(\\d+)_chunk.*", "\\1", basename(chunk_path)),
    ChunkID = sub(".*chunk(\\d+)", "\\1", basename(chunk_path)),
    GeneStart = NA_integer_,
    GeneEnd = NA_integer_,
    NumGenes = nrow(count_mat),
    NumPairs = nrow(count_mat) * (nrow(count_mat) - 1) / 2,
    GlobalPairs = global_pairs,
    BonferroniThreshold = bonf_thr,
    SignificantPairs = nrow(pairs_dt)
  )
  fwrite(summ_dt, file.path(chunk_path, "summary.tsv"), sep = "\t")
  nrow(pairs_dt)
}

tot_sig <- 0L
for (chunk in chunk_dirs) {
  n_sig <- process_chunk(chunk)
  if (!is.null(n_sig)) {
    tot_sig <- tot_sig + n_sig
  }
}

cat(sprintf("\nReapplied threshold for %d chunks. Total significant pairs: %d\n",
            length(chunk_dirs), tot_sig))
cat("Next steps:\n")
cat("  1. bash submit_step2b_merge.sh\n")
cat("  2. bash submit_step3.sh\n")
cat("  3. /apps/lib/R/R-4.4.0/bin/Rscript step4_merge_clusters.R\n\n")
