# MERFISH Louvain and aSBM real-data baselines using PC plus coordinate inputs.
# The script directly reuses the processed PCs and coordinates stored in
# coord_and_pc_19.RData and does not recompute PCA inside the run script.

script_file <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(dirname(script_file), "realdata_submission_helpers.R"))

paths <- submission_paths(script_file)
source("FitAttribute.R")

suppressPackageStartupMessages({
  library(FNN)
  library(igraph)
  library(mclust)
  library(mvtnorm)
})

data_obj <- load_merfish_processed_inputs(paths$data_root)
common_cells <- data_obj$cells
true_labels <- data_obj$true_labels

n_runs <- 10L
k_neighbors <- 6L

cat(sprintf("Loaded %d MERFISH cells from processed inputs.\n", length(common_cells)))

knn_res <- FNN::get.knn(data_obj$pca, k = k_neighbors)
edge_vector <- c()
for (i in seq_len(nrow(knn_res$nn.index))) {
  edge_vector <- c(edge_vector, as.vector(rbind(rep(i, k_neighbors), knn_res$nn.index[i, ])))
}
graph_obj <- simplify(make_graph(edge_vector, directed = FALSE))

louvain_metrics <- data.frame(Run = integer(), ARI = numeric(), NMI = numeric())
louvain_labels <- data.frame(Cell = common_cells, GroundTruth = true_labels, stringsAsFactors = FALSE)

for (run_id in seq_len(n_runs)) {
  set.seed(123 + run_id)
  louvain_out <- cluster_louvain(graph_obj)
  pred_labels <- membership(louvain_out)

  louvain_metrics <- rbind(
    louvain_metrics,
    data.frame(
      Run = run_id,
      ARI = adjustedRandIndex(true_labels, pred_labels),
      NMI = igraph::compare(as.numeric(factor(true_labels)), as.numeric(factor(pred_labels)), method = "nmi")
    )
  )
  louvain_labels[[paste0("Run", run_id)]] <- pred_labels
}

adj_matrix <- build_symmetric_knn_adjacency(data_obj$pca, k_neighbors)
attribute_data <- data_obj$spatial

asbm_metrics <- data.frame(Run = integer(), ARI = numeric(), NMI = numeric())
asbm_labels <- data.frame(Cell = common_cells, GroundTruth = true_labels, stringsAsFactors = FALSE)

for (run_id in seq_len(n_runs)) {
  set.seed(99 + run_id)
  fit_res <- FitAttribute(Network = adj_matrix, AttributeMat = attribute_data, Prob = 0)
  pred_labels <- fit_res$Comm

  asbm_metrics <- rbind(
    asbm_metrics,
    data.frame(
      Run = run_id,
      ARI = adjustedRandIndex(true_labels, pred_labels),
      NMI = igraph::compare(as.numeric(factor(true_labels)), as.numeric(factor(pred_labels)), method = "nmi")
    )
  )
  asbm_labels[[paste0("Run", run_id)]] <- pred_labels
}

write_realdata_outputs(paths$results_dir, "merfish_louvain_pc_coord_input_10runs", louvain_metrics, louvain_labels)
write_realdata_outputs(paths$results_dir, "merfish_asbm_pc_coord_input_10runs", asbm_metrics, asbm_labels)

cat("Finished MERFISH Louvain and aSBM processed-input runs.\n")
