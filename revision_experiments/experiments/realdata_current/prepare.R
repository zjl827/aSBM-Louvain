#!/usr/bin/env Rscript
# From the release root. Rebuild the fixed inputs from public processed PCs
# and the corresponding benchmark H5ADs.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(TRUE)
flag_indices <- 2L * seq_len(length(args) %/% 2L) - 1L
flags <- args[flag_indices]
if (length(args) %% 2L || anyDuplicated(flags) ||
    any(!flags %in% c("--input-dir", "--output-dir", "--processed-dir"))) {
  stop("Usage: Rscript experiments/realdata_current/prepare.R [--input-dir INPUT_DIR] [--output-dir DIR] [--processed-dir PC_DIR]")
}
options <- setNames(args[flag_indices + 1L], flags)
input <- if ("--input-dir" %in% names(options)) options[["--input-dir"]] else Sys.getenv("REVISION_H5AD_DIR", "data/realdata")
output <- if ("--output-dir" %in% names(options)) options[["--output-dir"]] else Sys.getenv("REVISION_REALDATA_DIR", "results/realdata_current/inputs")
processed <- if ("--processed-dir" %in% names(options)) options[["--processed-dir"]] else Sys.getenv("REVISION_PROCESSED_DIR", "..")
script_dir <- "experiments/realdata_current"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
tmp <- tempfile("revision-realdata-")
dir.create(tmp)

star_file <- file.path(processed, "starmap_PCA_scores_50PCs.csv")
pca_star <- as.matrix(read.csv(star_file, row.names = 1L, check.names = FALSE))
star_annotation <- read.delim(file.path(input, "Annotation_STARmap_20180505_BY3_1k.txt"),
  row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(setequal(rownames(star_annotation), rownames(pca_star)))
star_annotation <- star_annotation[rownames(pca_star), , drop = FALSE]
mer_env <- new.env(parent = emptyenv())
load(file.path(processed, "coord_and_pc_19.RData"), envir = mer_env)
pca_mer <- as.matrix(mer_env$gene_data_pc)
if (ncol(pca_mer) > nrow(pca_mer)) pca_mer <- t(pca_mer)
mer_file <- file.path(tmp, "merfish_pca_50.csv")
write.csv(pca_mer, mer_file, row.names = TRUE, quote = FALSE)
status <- system2(Sys.getenv("REVISION_PYTHON", "python3"),
  shQuote(c(file.path(script_dir, "prepare_neighbors.py"), input, star_file, mer_file, tmp, output)))
if (status != 0L) stop("Metadata and neighbor preparation failed.")

build_adjacency <- function(neighbors) {
  n <- nrow(neighbors)
  row_ids <- rep(seq_len(n), each = ncol(neighbors))
  column_ids <- as.integer(t(neighbors))
  result <- sparseMatrix(i = c(row_ids, column_ids), j = c(column_ids, row_ids), x = 1, dims = c(n, n))
  result@x[] <- 1
  diag(result) <- 0
  drop0(result)
}

for (dataset in c("STARmap", "MERFISH")) {
  stem <- tolower(dataset)
  pca <- if (dataset == "STARmap") pca_star else pca_mer
  metadata <- read.csv(file.path(tmp, paste0(stem, "_metadata.csv")), row.names = 1L,
    check.names = FALSE, stringsAsFactors = FALSE)
  cells <- rownames(pca)
  stopifnot(identical(rownames(metadata), cells), !anyNA(metadata$truth),
    nrow(pca) == if (dataset == "STARmap") 1207L else 5803L,
    ncol(pca) == 50L)
  spatial <- if (dataset == "STARmap") as.matrix(star_annotation[, c("X", "Y")]) else as.matrix(mer_env$coord)
  truth <- if (dataset == "STARmap") star_annotation$Annotation else metadata$truth
  stopifnot(length(truth) == nrow(pca), !anyNA(truth),
    length(unique(truth)) == if (dataset == "STARmap") 16L else 15L)
  stopifnot(nrow(spatial) == nrow(pca), ncol(spatial) == 2L, all(is.finite(spatial)))
  resolution <- if (dataset == "STARmap") 1.5 else 0.5
  for (neighbor_type in c("expression", "spatial")) {
    neighbors <- as.matrix(read.csv(file.path(tmp, paste0(stem, "_", neighbor_type, "_knn6.csv")),
      row.names = 1L, check.names = FALSE))
    stopifnot(identical(rownames(neighbors), cells), ncol(neighbors) == 6L)
    primary <- neighbor_type == "expression"
    orientation <- if (primary) "expression_graph_spatial_attributes" else "spatial_graph_expression_attributes"
    object <- list(dataset = dataset, orientation = orientation, adjacency = build_adjacency(neighbors),
      attributes = scale(if (primary) spatial else pca), truth = truth,
      resolution = resolution, graph_definition = paste0(neighbor_type,
        if (primary) "_pca_symmetric_knn6" else "_symmetric_knn6"),
      attribute_definition = if (primary) "scaled_spatial_coordinates" else "scaled_expression_pca_50",
      cell_ids = cells)
    path <- file.path(output, paste0(dataset, "__", orientation, ".rds"))
    saveRDS(object, path)
    message("Wrote ", path)
  }
}
unlink(tmp, recursive = TRUE)
