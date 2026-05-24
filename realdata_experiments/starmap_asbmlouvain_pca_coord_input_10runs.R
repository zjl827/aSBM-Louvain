# STARmap aSBM-Louvain real-data experiment using PCA-score plus coordinate
# inputs. The script directly reuses the precomputed STARmap PCA file and does
# not recompute PCA inside the run script.

script_file <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(dirname(script_file), "realdata_submission_helpers.R"))

paths <- submission_paths(script_file)
source("asbmlou.R")

suppressPackageStartupMessages({
  library(FNN)
  library(igraph)
  library(mclust)
  library(mvtnorm)
})

data_obj <- load_starmap_processed_inputs(paths$data_root)
common_cells <- data_obj$cells
true_labels <- data_obj$true_labels
coords_mat <- data_obj$spatial

n_runs <- 10L
k_neighbors <- 6L
target_resolution <- 1.5

cat(sprintf("Loaded %d STARmap cells from processed inputs.\n", length(common_cells)))

adj_matrix <- build_symmetric_knn_adjacency(data_obj$pca, k_neighbors)
attribute_data <- scale(coords_mat)

metrics_df <- data.frame(
  Run = integer(),
  K_Found = integer(),
  ARI = numeric(),
  NMI = numeric()
)
labels_df <- data.frame(Cell = common_cells, GroundTruth = true_labels, stringsAsFactors = FALSE)

for (run_id in seq_len(n_runs)) {
  set.seed(42 + run_id)
  fit_res <- aSBM_Louvain(
    AdjMat = adj_matrix,
    AttributeMat = attribute_data,
    rho_tau = 0.8,
    T_outer = 1000,
    init_method = "kmeans",
    resolution = target_resolution
  )
  pred_labels <- fit_res$partition

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

write_realdata_outputs(paths$results_dir, "starmap_asbmlouvain_pca_coord_input_10runs", metrics_df, labels_df)

cat("Finished STARmap aSBM-Louvain processed-input runs.\n")
