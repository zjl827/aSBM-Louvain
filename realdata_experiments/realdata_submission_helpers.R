# Shared helpers for the clean real-data submission scripts.
# These helpers preserve the original input definition used by each method while
# giving the cleaned scripts a consistent loading and output interface.

this_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    stop("Unable to determine the current script path from commandArgs().")
  }
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
}

submission_paths <- function(current_file = this_script_path()) {
  script_file <- normalizePath(current_file, winslash = "/", mustWork = TRUE)
  submission_dir <- dirname(script_file)
  project_root <- normalizePath(file.path(submission_dir, ".."), winslash = "/", mustWork = TRUE)
  data_root <- project_root
  results_dir <- file.path(submission_dir, "results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  list(
    script_file = script_file,
    submission_dir = submission_dir,
    project_root = project_root,
    data_root = data_root,
    results_dir = results_dir
  )
}

load_starmap_processed_inputs <- function(data_root) {
  pca_scores <- read.csv(file.path(data_root, "starmap_PCA_scores_50PCs.csv"), row.names = 1)
  spatial <- read.csv(file.path(data_root, "starmap_spatial.csv"), row.names = 1)
  annotation <- read.table(
    file.path(data_root, "Annotation_STARmap_20180505_BY3_1k.txt"),
    header = TRUE,
    sep = "\t",
    row.names = 1
  )

  common_cells <- Reduce(intersect, list(rownames(pca_scores), rownames(spatial), rownames(annotation)))
  common_cells <- sort(common_cells)
  if (length(common_cells) == 0) {
    stop("No overlapping STARmap cell IDs were found across PCA/spatial/annotation.")
  }

  spatial_sub <- spatial[common_cells, c("X", "Y"), drop = FALSE]
  list(
    cells = common_cells,
    pca = as.matrix(pca_scores[common_cells, , drop = FALSE]),
    spatial = as.matrix(spatial_sub),
    true_labels = annotation[common_cells, "Annotation", drop = TRUE]
  )
}

load_merfish_processed_inputs <- function(data_root) {
  input_env <- new.env(parent = emptyenv())
  load(file.path(data_root, "coord_and_pc_19.RData"), envir = input_env)
  annotation <- read.csv(
    file.path(data_root, "cell_type_annotation_0.19.csv"),
    row.names = 1,
    check.names = FALSE
  )

  if (ncol(input_env$gene_data_pc) > nrow(input_env$gene_data_pc)) {
    pca_cells_by_pc <- t(input_env$gene_data_pc)
  } else {
    pca_cells_by_pc <- input_env$gene_data_pc
  }

  if (nrow(input_env$coord) != nrow(pca_cells_by_pc)) {
    stop("MERFISH PCA and coordinate objects have inconsistent row counts.")
  }
  rownames(input_env$coord) <- rownames(pca_cells_by_pc)

  common_cells <- Reduce(intersect, list(rownames(pca_cells_by_pc), rownames(input_env$coord), rownames(annotation)))
  common_cells <- sort(common_cells)
  if (length(common_cells) == 0) {
    stop("No overlapping MERFISH cell IDs were found across processed inputs and annotation.")
  }

  coord_sub <- input_env$coord[common_cells, , drop = FALSE]
  if (!all(c("coord_x", "coord_y") %in% colnames(coord_sub))) {
    colnames(coord_sub)[1:2] <- c("coord_x", "coord_y")
  }

  bact_gene_data_pc <- if (ncol(input_env$gene_data_pc) > nrow(input_env$gene_data_pc)) {
    input_env$gene_data_pc[, common_cells, drop = FALSE]
  } else {
    t(input_env$gene_data_pc[common_cells, , drop = FALSE])
  }

  list(
    cells = common_cells,
    pca = as.matrix(pca_cells_by_pc[common_cells, , drop = FALSE]),
    spatial = as.matrix(coord_sub[, c("coord_x", "coord_y"), drop = FALSE]),
    true_labels = annotation[common_cells, 1, drop = TRUE],
    bact_gene_data_pc = as.matrix(bact_gene_data_pc)
  )
}

build_symmetric_knn_adjacency <- function(feature_mat, k) {
  knn_res <- FNN::get.knn(as.matrix(feature_mat), k = k)
  n_obs <- nrow(knn_res$nn.index)
  adj_matrix <- matrix(0, nrow = n_obs, ncol = n_obs)
  for (i in seq_len(n_obs)) {
    neighbor_idx <- knn_res$nn.index[i, ]
    adj_matrix[i, neighbor_idx] <- 1
    adj_matrix[neighbor_idx, i] <- 1
  }
  adj_matrix
}

get_mode <- function(values) {
  uniq_values <- unique(values)
  uniq_values[which.max(tabulate(match(values, uniq_values)))]
}

metric_summary_frame <- function(metrics_df) {
  numeric_cols <- setdiff(names(metrics_df)[vapply(metrics_df, is.numeric, logical(1))], "Run")
  summary_list <- list(runs = nrow(metrics_df))
  for (col_name in numeric_cols) {
    col_values <- metrics_df[[col_name]]
    summary_list[[paste0("mean_", col_name)]] <- mean(col_values, na.rm = TRUE)
    summary_list[[paste0("sd_", col_name)]] <- sd(col_values, na.rm = TRUE)
    summary_list[[paste0("median_", col_name)]] <- median(col_values, na.rm = TRUE)
  }
  as.data.frame(summary_list, check.names = FALSE)
}

write_realdata_outputs <- function(results_dir, prefix, metrics_df, labels_df) {
  metrics_path <- file.path(results_dir, paste0(prefix, "_metrics.csv"))
  labels_path <- file.path(results_dir, paste0(prefix, "_labels.csv"))
  summary_path <- file.path(results_dir, paste0(prefix, "_summary.csv"))

  write.csv(metrics_df, metrics_path, row.names = FALSE)
  write.csv(labels_df, labels_path, row.names = FALSE)
  write.csv(metric_summary_frame(metrics_df), summary_path, row.names = FALSE)

  invisible(list(
    metrics = metrics_path,
    labels = labels_path,
    summary = summary_path
  ))
}
