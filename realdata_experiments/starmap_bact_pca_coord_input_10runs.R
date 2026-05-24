# STARmap BACT real-data experiment using PCA-score plus coordinate inputs.
# The BACT input gene_data_pc is built directly from the precomputed STARmap
# PCA file, so the script does not recompute PCA inside the run script.

script_file <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(dirname(script_file), "realdata_submission_helpers.R"))

paths <- submission_paths(script_file)

suppressPackageStartupMessages({
  library(BACT)
  library(igraph)
  library(mclust)
})

data_obj <- load_starmap_processed_inputs(paths$data_root)
common_cells <- data_obj$cells
true_labels <- data_obj$true_labels

n_runs <- 10L

coord_use <- as.data.frame(data_obj$spatial)
colnames(coord_use) <- c("coord_x", "coord_y")
gene_data_pc <- t(data_obj$pca)

metrics_df <- data.frame(Run = integer(), K_Found = integer(), ARI = numeric(), NMI = numeric())
labels_df <- data.frame(Cell = common_cells, GroundTruth = true_labels, stringsAsFactors = FALSE)

for (run_id in seq_len(n_runs)) {
  seed_i <- 99 + run_id
  fit_res <- BACT(
    gene_data_pc = gene_data_pc,
    coord = coord_use,
    platform = "sc",
    num_init = 7,
    num_nei = 6,
    a_eta = 0,
    b_eta = 1.5,
    IGkappa = 2,
    IGtau = 10,
    dpAlpha = 1,
    a_beta = 1,
    tau_beta = 1,
    tau0 = 0.01,
    tau1 = 0.05,
    M0 = 50,
    numOfMCMC = 6000,
    burnIn = 3000,
    Is_beta_zero = FALSE,
    Is_warm_start = TRUE,
    Is_kmeans_use_mean_sd = TRUE,
    Is_print = TRUE,
    print_gap = 500,
    Is_random_seed = TRUE,
    random_seed = seed_i
  )

  cl_ids_mcmc <- fit_res$clIds_mcmc
  pred_labels <- if (ncol(cl_ids_mcmc) == length(common_cells)) {
    apply(cl_ids_mcmc, 2, get_mode)
  } else if (nrow(cl_ids_mcmc) == length(common_cells)) {
    apply(cl_ids_mcmc, 1, get_mode)
  } else {
    stop("clIds_mcmc dimensions do not match the STARmap cell count.")
  }

  metrics_df <- rbind(
    metrics_df,
    data.frame(
      Run = run_id,
      K_Found = length(unique(pred_labels)),
      ARI = adjustedRandIndex(true_labels, pred_labels),
      NMI = igraph::compare(as.numeric(factor(true_labels)), as.numeric(factor(pred_labels)), method = "nmi")
    )
  )
  labels_df[[paste0("Run", run_id)]] <- pred_labels
}

write_realdata_outputs(paths$results_dir, "starmap_bact_pca_coord_input_10runs", metrics_df, labels_df)

cat("Finished STARmap BACT processed-input runs.\n")
