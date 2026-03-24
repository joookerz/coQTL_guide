# coQTL guide

This directory contains a simplified coQTL workflow for one cell type. A user provides one single-cell count matrix, one gene annotation file, and genotype files for SAIGE-QTL, then runs the full pipeline step by step to generate one complete result set.

## What this workflow does

1. identify co-expression clusters from a single-cell count matrix
2. perform PCA for each inferred cluster
3. prepare cluster-level phenotype files for SAIGE-QTL
4. run SAIGE-QTL to generate association results

## Files in this directory

- `bin/`: command-line entry points
- `templates/config.template.R`: config template users edit for their own data
- `templates/example_config.R`: example config with filled path patterns
- `scripts/method2/`: cluster-identification scripts
- `scripts/pcqtl/`: PCA and SAIGE input-preparation scripts
- `envs/coqtl_r.yml`: base R environment for Method 2 and PCA

## Software requirements

### Method 2 and PCA

```bash
conda env create -f envs/coqtl_r.yml
conda activate coqtl-r
Rscript -e "install.packages(c('data.table','optparse','remotes'), repos='https://cloud.r-project.org')"
```

Required R packages:

- `data.table`
- `optparse`
- `fasthurdle >= 1.1.1`

`envs/coqtl_r.yml` already covers the base R environment used by this guide. `fasthurdle` is the only required R package that still needs to be installed separately before running Method 2.

Example checks:

```bash
Rscript -e "packageVersion('data.table'); packageVersion('optparse')"
Rscript -e "packageVersion('fasthurdle')"
```

### SAIGE-QTL

This guide does not include SAIGE-QTL installation instructions.

Please use the official SAIGE-QTL documentation:

- https://weizhou0.github.io/SAIGE-QTL-doc/

After SAIGE-QTL is installed, set the related command fields in `config.R`, such as:

- `SAIGE_CMD_PREFIX`
- `SAIGE_STEP1_CMD`
- `SAIGE_STEP2_CMD`
- `SAIGE_STEP3_CMD`

## Required input files

### 1. Single-cell count matrix

Set as `COUNT_FILE` in `config.R`.

Expected format:

- tab-delimited `.tsv` or `.tsv.gz`
- one row per cell
- one column per gene
- must contain a cell ID column: `CellID` or `barcode`
- must contain an individual/sample column: `individual` or `IndividualID`
- may contain covariate columns such as `age`, `sex`, `pc1`-`pc6`, `pf1`, `pf2`

All remaining non-metadata columns are treated as gene-expression columns.

### 2. Gene annotation file

Set as `GENE_INFO_FILE` in `config.R`.

Required columns:

- `gene_name`
- `chr_numeric`
- `start`
- `end`

### 3. Genotype input for SAIGE-QTL

Required only for the SAIGE stage.

Set these fields in `config.R`:

- `VR_PLINK_PREFIX`
- `GENO_BED_PATTERN`
- `GENO_BIM_PATTERN`
- `GENO_FAM_PATTERN`

Use `{chr}` as the chromosome placeholder in the genotype filename patterns.

## Minimal usage

### Step 1. Create a work directory

```bash
bash bin/setup_celltype.sh B_in /path/to/coqtl_runs/B_in
```

This creates a new work directory containing:

- `config.R`
- `cluster_identification/method2_sc_hurdle/`
- `pcQTL/`

### Step 2. Edit `config.R`

The user must fill in at least:

- `CELL_TYPE`
- `COUNT_FILE`
- `GENE_INFO_FILE`
- `WORK_ROOT`
- `HURDLE_COVAR_COLS` if hurdle models should adjust for covariates
- `COVAR_COLS`
- `SAIGE_COVAR_COLS`
- `VR_PLINK_PREFIX`
- `GENO_BED_PATTERN`
- `GENO_BIM_PATTERN`
- `GENO_FAM_PATTERN`
- `SAIGE_CMD_PREFIX` if SAIGE-QTL is run through Docker or Singularity

An example filled template is available at `templates/example_config.R`.

To add or remove covariates, modify:

- `HURDLE_COVAR_COLS`
- `COVAR_COLS`
- `SAIGE_COVAR_COLS`

The corresponding columns must exist in the count matrix.

### Hurdle covariates

Method 2 step 2 supports optional covariate adjustment through `HURDLE_COVAR_COLS` in `config.R`.

- if `HURDLE_COVAR_COLS <- character(0)`, the hurdle model uses no extra covariates
- if `HURDLE_COVAR_COLS <- c("age","sex","pc1","pc2")`, each pairwise hurdle fit adjusts for those columns from `COUNT_FILE`

Current model form is:

- `count_i ~ count_j + covariates`
- `count_j ~ count_i + covariates`

The count matrix must contain every column listed in `HURDLE_COVAR_COLS`.

### Step 3. Run the complete workflow

```bash
bash bin/run_method2.sh /path/to/coqtl_runs/B_in/config.R
bash bin/run_pcqtl_pca.sh /path/to/coqtl_runs/B_in/config.R
bash bin/run_saige_all.sh /path/to/coqtl_runs/B_in/config.R
```

These three commands are the main user workflow.

## What each command produces

### `run_method2.sh`

This step:

1. filters sparse genes
2. computes gene-gene hurdle associations in chunks
3. merges significant pairs
4. identifies co-expression clusters
5. merges cluster results across chromosomes

Main outputs:

- `cluster_identification/results/method2_sc_hurdle/filtered_genes.tsv`
- `cluster_identification/results/method2_sc_hurdle/gene_associations/`
- `cluster_identification/results/method2_sc_hurdle/clusters_by_chr/`
- `cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_summary.tsv`
- `cluster_identification/results/method2_sc_hurdle/merged_clusters/cluster_gene_assignments.tsv`

### `run_pcqtl_pca.sh`

This step:

1. performs PCA for each inferred cluster
2. writes phenotype files containing covariates and cluster PCs
3. prepares SAIGE region files and the cluster-PC map

Main outputs:

- `pcQTL/step2_pca/<cluster_id>/pheno_with_pcs.tsv`
- `pcQTL/step2_pca/<cluster_id>/variance_explained.tsv`
- `pcQTL/step2_pca/<cluster_id>/gene_loadings.tsv`
- `pcQTL/step2_pca/all_clusters_summary.tsv`
- `pcQTL/step2_pca/pc_list_for_saige.tsv`
- `pcQTL/step3_saige/regions/`
- `pcQTL/step3_saige/cluster_pc_map.tsv`

### `run_saige_all.sh`

This step runs SAIGE-QTL sequentially for all cluster-PC phenotypes prepared in the previous step.

Main outputs:

- `pcQTL/step3_saige/step1/`
- `pcQTL/step3_saige/step2/`
- `pcQTL/step3_saige/step3/`

If the user wants to test only one cluster first:

```bash
bash bin/run_saige_cluster.sh /path/to/coqtl_runs/B_in/config.R SC_chr10_cluster_001
```

### SAIGE-QTL defaults used in this guide

The wrapper `bin/run_saige_cluster.sh` currently uses these default SAIGE-QTL settings:

- `--sampleIDColinphenoFile=individual`
- `--traitType=quantitative`
- `--invNormalize=TRUE`
- `--useSparseGRMtoFitNULL=FALSE`
- `--useGRMtoFitNULL=FALSE`
- `--skipVarianceRatioEstimation=FALSE`
- `--isRemoveZerosinPheno=FALSE`
- `--isCovariateOffset=FALSE`
- `--isCovariateTransform=TRUE`
- `--skipModelFitting=FALSE`
- `--tol=0.00001`
- `--LOCO=FALSE`
- `--minMAF=0.05`
- `--markers_per_chunk=10000`
- cis-window for each cluster region: `500000` bp on each side, controlled by `CIS_WINDOW`

Required phenotype/covariate columns used by the wrapper:

- `individual`
- all columns listed in `SAIGE_COVAR_COLS`
- cluster PC columns generated in `pheno_with_pcs.tsv`, such as `PC1`, `PC2`, ...

In the example configuration used here, the SAIGE covariates are:

- `age`
- `sex`
- `pc1`
- `pc2`
- `pc3`
- `pc4`
- `pc5`
- `pc6`
- `pf1`
- `pf2`

## Final result directory

All outputs for one cell type are written under `WORK_ROOT` in `config.R`.

The final result set includes:

- inferred co-expression clusters
- cluster-level PCA phenotypes
- SAIGE null models
- SAIGE variant-level association results
- SAIGE gene-level p-value results

## Notes

- The guide is intentionally sequential and focused on producing one complete result.
- It does not include Slurm instructions.
- The Method 2 and PCA environment is complete except for `fasthurdle`, which must be installed by the user.
- SAIGE-QTL installation and runtime setup should follow the official documentation: https://weizhou0.github.io/SAIGE-QTL-doc/
