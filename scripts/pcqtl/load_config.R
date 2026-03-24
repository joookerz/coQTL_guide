get_script_path <- function() {
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

load_pcqtl_config <- function(script_dir = NULL) {
  config_path <- find_config_path(script_dir)
  cfg_env <- new.env()
  sys.source(config_path, envir = cfg_env)

  required <- c(
    "CELL_TYPE", "COUNT_FILE", "GENE_INFO_FILE", "WORK_ROOT",
    "COVAR_COLS", "VARIANCE_THRESHOLD", "CIS_WINDOW",
    "SAIGE_COVAR_COLS", "SAIGE_STEP1_CMD", "SAIGE_STEP2_CMD",
    "SAIGE_STEP3_CMD", "SAIGE_CMD_PREFIX", "VR_PLINK_PREFIX",
    "GENO_BED_PATTERN", "GENO_BIM_PATTERN", "GENO_FAM_PATTERN",
    "MIN_MAF", "MARKERS_PER_CHUNK"
  )
  missing <- required[!vapply(required, exists, logical(1), envir = cfg_env)]
  if (length(missing)) {
    stop("Config file missing definitions: ", paste(missing, collapse = ", "))
  }

  cell_id_candidates <- if (exists("CELL_ID_COL_CANDIDATES", envir = cfg_env)) {
    cfg_env$CELL_ID_COL_CANDIDATES
  } else {
    c("CellID", "barcode")
  }
  individual_candidates <- if (exists("INDIVIDUAL_COL_CANDIDATES", envir = cfg_env)) {
    cfg_env$INDIVIDUAL_COL_CANDIDATES
  } else {
    c("individual", "IndividualID")
  }
  celltype_candidates <- if (exists("CELLTYPE_COL_CANDIDATES", envir = cfg_env)) {
    cfg_env$CELLTYPE_COL_CANDIDATES
  } else {
    c("CellType")
  }

  pcqtl_dir <- file.path(cfg_env$WORK_ROOT, "pcQTL")
  cluster_file <- file.path(
    cfg_env$WORK_ROOT,
    "cluster_identification", "results", "method2_sc_hurdle",
    "merged_clusters", "cluster_summary.tsv"
  )
  cluster_gene_file <- file.path(
    cfg_env$WORK_ROOT,
    "cluster_identification", "results", "method2_sc_hurdle",
    "merged_clusters", "cluster_gene_assignments.tsv"
  )

  list(
    config_path = config_path,
    cell_type = cfg_env$CELL_TYPE,
    count_file = cfg_env$COUNT_FILE,
    gene_info_file = cfg_env$GENE_INFO_FILE,
    work_root = cfg_env$WORK_ROOT,
    pcqtl_dir = pcqtl_dir,
    cluster_file = cluster_file,
    cluster_gene_file = cluster_gene_file,
    covar_cols = cfg_env$COVAR_COLS,
    variance_threshold = cfg_env$VARIANCE_THRESHOLD,
    cis_window = cfg_env$CIS_WINDOW,
    saige_covar_cols = cfg_env$SAIGE_COVAR_COLS,
    saige_step1_cmd = cfg_env$SAIGE_STEP1_CMD,
    saige_step2_cmd = cfg_env$SAIGE_STEP2_CMD,
    saige_step3_cmd = cfg_env$SAIGE_STEP3_CMD,
    saige_cmd_prefix = cfg_env$SAIGE_CMD_PREFIX,
    vr_plink_prefix = cfg_env$VR_PLINK_PREFIX,
    geno_bed_pattern = cfg_env$GENO_BED_PATTERN,
    geno_bim_pattern = cfg_env$GENO_BIM_PATTERN,
    geno_fam_pattern = cfg_env$GENO_FAM_PATTERN,
    min_maf = cfg_env$MIN_MAF,
    markers_per_chunk = cfg_env$MARKERS_PER_CHUNK,
    cell_id_candidates = cell_id_candidates,
    individual_col_candidates = individual_candidates,
    celltype_col_candidates = celltype_candidates
  )
}
