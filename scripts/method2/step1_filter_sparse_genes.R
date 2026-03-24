#!/usr/bin/env Rscript

# =============================================================================
# Method 2 - Step 1: Filter sparse genes using 1% nonzero cutoff
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
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

cat(sprintf("\n=== Method 2 Step 1: Filter sparse genes (%s) ===\n", cfg$cell_type))
cat(sprintf("Start time: %s\n\n", Sys.time()))

COUNT_FILE <- cfg$count_file
GENE_INFO_FILE <- cfg$gene_info_file
OUTPUT_DIR <- cfg$method2_results
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

NONZERO_CUTOFF <- cfg$nonzero_cutoff

cat("Loading count matrix (this may take a while)...\n")
count_dt <- fread(COUNT_FILE)
id_col <- intersect(c("CellID", "barcode"), names(count_dt))
if (!length(id_col)) stop("Count file must contain CellID or barcode column")
id_col <- id_col[1]
genes <- setdiff(names(count_dt), c(id_col, "IndividualID", "individual", "CellType"))
n_cells <- nrow(count_dt)
cat(sprintf("  Cells: %d  Genes: %d\n", n_cells, length(genes)))

cat("Calculating nonzero proportions...\n")
count_mat <- as.matrix(count_dt[, ..genes])
nonzero_prop <- colMeans(count_mat > 0)
rm(count_mat); gc()

nonzero_dt <- data.table(
  gene_name = genes,
  nonzero_prop = nonzero_prop
)

gene_info <- fread(GENE_INFO_FILE, select = c("gene_name", "chr_numeric", "start", "end"))
result_dt <- merge(nonzero_dt, gene_info, by = "gene_name", all.x = TRUE)
result_dt[, keep := nonzero_prop >= NONZERO_CUTOFF]

cat(sprintf("  Genes passing cutoff: %d / %d (%.1f%%)\n",
            sum(result_dt$keep), nrow(result_dt),
            100 * mean(result_dt$keep)))

out_file <- file.path(OUTPUT_DIR, "filtered_genes.tsv")
fwrite(result_dt, out_file, sep = "\t")

cat("\n=== Complete ===\n")
cat(sprintf("Saved: %s\n", out_file))
cat(sprintf("End time: %s\n", Sys.time()))
