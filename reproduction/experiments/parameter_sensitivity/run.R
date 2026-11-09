#!/usr/bin/env Rscript
# Supplement S11: one-factor-at-a-time sweeps on the fixed real datasets.
suppressPackageStartupMessages({
  library(igraph)
  library(Matrix)
  library(mclust)
  library(mvtnorm)
})

script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
suite <- dirname(script)
code <- normalizePath(file.path(suite, "../../code"))
source(file.path(code, "asbmlou.R"))
source(file.path(code, "fit_asbm_louvain.R"))
source(file.path(code, "pure_asbm_vectorized.R"))

args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript experiments/parameter_sensitivity/run.R OUTPUT_DIR [CASE_ID ...]")
out <- args[1L]
real_dir <- Sys.getenv("REVISION_REALDATA_DIR", "results/realdata_current/inputs")
knn_dir <- Sys.getenv("REVISION_KNN_DIR", real_dir)
tasks <- read.csv(file.path(suite, "tasks.csv"))
stopifnot(nrow(tasks) == 340L, !anyDuplicated(tasks$key))
if (length(args) > 1L) {
  stopifnot(all(args[-1L] %in% tasks$key))
  tasks <- tasks[tasks$key %in% args[-1L], ]
}

load_dataset <- function(dataset) {
  input_path <- file.path(real_dir, paste0(dataset, "__expression_graph_spatial_attributes.rds"))
  neighbor_path <- file.path(knn_dir, paste0(tolower(dataset), "_expression_knn12.csv"))
  if (!all(file.exists(input_path, neighbor_path))) {
    stop("Missing expression-graph inputs for ", dataset, ". Run experiments/realdata_current/prepare.R.")
  }
  input <- readRDS(input_path)
  neighbors <- as.matrix(read.csv(neighbor_path, row.names = 1, check.names = FALSE))
  n <- nrow(input$attributes)
  stopifnot(nrow(neighbors) == n, ncol(neighbors) == 12L, length(input$truth) == n,
            all(dim(input$adjacency) == n), all(is.finite(neighbors)),
            all(neighbors == as.integer(neighbors)), all(neighbors >= 1L), all(neighbors <= n),
            all(neighbors != seq_len(n)), all(apply(neighbors, 1L, anyDuplicated) == 0L))
  cell_ids <- input$cell_ids
  if (is.null(cell_ids)) cell_ids <- rownames(input$attributes)
  if (!is.null(cell_ids)) stopifnot(identical(rownames(neighbors), as.character(cell_ids)))
  build_graph <- function(k) {
    A <- sparseMatrix(i = rep(seq_len(n), each = k),
                      j = as.integer(t(neighbors[, seq_len(k), drop = FALSE])),
                      x = 1, dims = c(n, n))
    A <- as((A + t(A)) > 0, "dMatrix")
    diag(A) <- 0
    dimnames(A) <- dimnames(input$adjacency)
    A
  }
  # Use the 6NN baseline to check alignment, including older caches without cell IDs.
  stopifnot(all(build_graph(6L) == input$adjacency))
  counts <- unique(tasks$k[tasks$dataset == dataset])
  input$graphs <- setNames(lapply(counts, build_graph), as.character(counts))
  input
}

datasets <- unique(tasks$dataset)
inputs <- setNames(lapply(datasets, load_dataset), datasets)
dir.create(file.path(out, "runs"), recursive = TRUE, showWarnings = FALSE)

run_case <- function(task) {
  d <- inputs[[task$dataset]]
  A <- d$graphs[[as.character(task$k)]]
  set.seed(task$seed)
  init <- pure_asbm_matched_kmeans_initialization(A, d$attributes, resolution = task$resolution, nstart = 5L)
  fit <- aSBM_Louvain_revised(
    A, d$attributes, T_outer = 1000L, rho_tau = task$rho_tau, tau0 = task$tau0,
    resolution = task$resolution, epsilon_Q = task$epsilon_Q, epsilon_L = task$epsilon_L,
    k_stable_required = task$k_stable_required, init_method = "kmeans",
    initial_partition = init$partition, prior_feedback = TRUE
  )
  data.frame(task, ARI = adjustedRandIndex(d$truth, fit$partition),
             NMI = compare(as.integer(factor(d$truth)), as.integer(factor(fit$partition)), method = "nmi"),
             initial_K = init$estimated_groups, final_K = fit$final_K, iterations = fit$iterations_run,
             stop_reason = fit$stop_reason, converged = fit$converged)
}

rows <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i, , drop = FALSE]
  cat(sprintf("[%d/%d] %s\n", i, nrow(tasks), task$key))
  metrics <- run_case(task)
  write.csv(metrics, file.path(out, "runs", paste0(task$key, ".csv")), row.names = FALSE)
  metrics
})
raw <- do.call(rbind, rows)
write.csv(raw, file.path(out, "raw_metrics.csv"), row.names = FALSE)
groups <- split(raw, list(raw$dataset, raw$config_id), drop = TRUE)
summary <- lapply(groups, function(x) {
  columns <- c("ARI", "NMI", "initial_K", "final_K", "iterations")
  averages <- setNames(lapply(x[columns], mean), paste0(columns, "_mean"))
  deviations <- setNames(lapply(x[columns], sd), paste0(columns, "_sd"))
  data.frame(x[1L, c("dataset", "family", "config_id", "parameter_label", "k", "resolution_multiplier",
                     "resolution", "rho_tau", "tau0", "epsilon_L", "epsilon_Q")],
             n = nrow(x), averages, deviations)
})
write.csv(do.call(rbind, summary), file.path(out, "configuration_summary.csv"), row.names = FALSE)
