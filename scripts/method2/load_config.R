get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_idx <- grep("--file=", cmd_args)
  if (length(file_idx) > 0) {
    return(normalizePath(sub("--file=", "", cmd_args[file_idx])))
  }
  # fallback for interactive sessions
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
  normalizePath(".")
}

find_config_path <- function(script_dir = NULL) {
  if (is.null(script_dir)) {
    script_dir <- dirname(get_script_path())
  }
  env_cfg <- Sys.getenv("COQTL_CONFIG", unset = "")
  candidates <- unique(c(
    env_cfg,
    file.path(script_dir, "config.R"),
    file.path(script_dir, "..", "config.R"),
    file.path(script_dir, "..", "..", "config.R")
  ))
  for (cfg in candidates) {
    if (nzchar(cfg) && file.exists(cfg)) {
      return(normalizePath(cfg))
    }
  }
  stop("Could not find config.R. Set COQTL_CONFIG or place config.R in the cell type root.")
}

load_config <- function(script_dir = NULL) {
  config_path <- find_config_path(script_dir)

  cfg_env <- new.env()
  sys.source(config_path, envir = cfg_env)

  required <- c(
    "CELL_TYPE", "COUNT_FILE", "GENE_INFO_FILE", "WORK_ROOT",
    "NONZERO_CUTOFF", "P_THRESHOLD", "CHUNK_SIZE_GENES", "SECS_PER_PAIR",
    "MAX_WINDOW_SIZE", "MIN_WINDOW_SIZE", "CLUSTER_THRESHOLD"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = cfg_env)]
  if (length(missing)) {
    stop("Config file missing definitions: ", paste(missing, collapse = ", "))
  }

  work_dir <- file.path(cfg_env$WORK_ROOT, "cluster_identification")
  results_root <- file.path(work_dir, "results")
  method2_results <- file.path(results_root, "method2_sc_hurdle")
  logs_dir <- file.path(work_dir, "method2_sc_hurdle", "logs")
  dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(method2_results, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  hurdle_covar_cols <- if (exists("HURDLE_COVAR_COLS", envir = cfg_env)) {
    cfg_env$HURDLE_COVAR_COLS
  } else {
    character(0)
  }

  cfg <- list(
    config_path = config_path,
    cell_type = cfg_env$CELL_TYPE,
    work_root = cfg_env$WORK_ROOT,
    work_dir = work_dir,
    count_file = cfg_env$COUNT_FILE,
    gene_info_file = cfg_env$GENE_INFO_FILE,
    results_root = results_root,
    method2_results = method2_results,
    logs_dir = logs_dir,
    nonzero_cutoff = cfg_env$NONZERO_CUTOFF,
    p_threshold = cfg_env$P_THRESHOLD,
    chunk_size_genes = cfg_env$CHUNK_SIZE_GENES,
    secs_per_pair = cfg_env$SECS_PER_PAIR,
    max_window_size = cfg_env$MAX_WINDOW_SIZE,
    min_window_size = cfg_env$MIN_WINDOW_SIZE,
    cluster_threshold = cfg_env$CLUSTER_THRESHOLD,
    hurdle_covar_cols = hurdle_covar_cols
  )
  cfg$chunk_dir <- file.path(cfg$method2_results, "gene_associations_chunked")
  cfg$assoc_dir <- file.path(cfg$method2_results, "gene_associations")
  cfg$clusters_dir <- file.path(cfg$method2_results, "clusters_by_chr")
  cfg
}
