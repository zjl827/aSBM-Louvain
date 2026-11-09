#!/usr/bin/env Rscript
# Supplement S2: 3:1 and 5:1 unequal-size synthetic community experiments.
suppressPackageStartupMessages({
  library(igraph)
  library(mclust)
  library(mvtnorm)
})

script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
suite <- dirname(script)
code <- normalizePath(file.path(suite, "../../code"))
model <- new.env(parent = .GlobalEnv)
sys.source(file.path(code, "asbmlou.R"), envir = model)
sys.source(file.path(code, "asbmlou_revised_stop.R"), envir = model)
source(file.path(code, "pure_asbm_vectorized.R"))
source(file.path(code, "simulation_design.R"))

args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript experiments/unequal_sizes/run.R OUTPUT_DIR [CASE_ID ...]")
out <- args[1L]
tasks <- read.csv(file.path(suite, "tasks.csv"))
stopifnot(nrow(tasks) == 360L, !anyDuplicated(tasks$case_id))
if (length(args) > 1L) {
  stopifnot(all(args[-1L] %in% tasks$case_id))
  tasks <- tasks[tasks$case_id %in% args[-1L], ]
}
dir.create(file.path(out, "runs"), recursive = TRUE, showWarnings = FALSE)

score_fit <- function(task, method, truth, partition, iterations = NA_integer_, reason = "reference_partition") {
  stopifnot(length(truth) == length(partition), !anyNA(partition))
  data.frame(task, method = method, ARI = adjustedRandIndex(truth, partition),
             NMI = compare(as.integer(factor(truth)), as.integer(factor(partition)), method = "nmi"),
             final_K = length(unique(partition)), iterations = iterations, stop_reason = reason)
}

run_case <- function(task) {
  sizes <- as.integer(strsplit(task$block_sizes, ":", fixed = TRUE)[[1L]])
  d <- simulate_unbalanced_data(sizes, task$D, task$r, task$sigma, task$target_degree, task$data_seed)
  set.seed(task$algorithm_seed)
  initial <- pure_asbm_matched_kmeans_initialization(d$A, d$X, resolution = 1, nstart = 5L)$partition
  post_init_rng <- .Random.seed
  assign(".Random.seed", post_init_rng, envir = .GlobalEnv)
  full <- model$aSBM_Louvain_revised(
    d$A, d$X, T_outer = 1000L, rho_tau = 0.8, tau0 = 1, resolution = 1,
    epsilon_Q = 1e-4, epsilon_L = 1e-3, k_stable_required = 3L,
    init_method = "kmeans", initial_partition = initial, prior_feedback = TRUE
  )
  assign(".Random.seed", post_init_rng, envir = .GlobalEnv)
  pure <- fit_pure_asbm_vectorized(
    d$A, d$X, initial_partition = initial, tolerance = 1, timeout_sec = 3600,
    max_iterations = 1000L
  )
  set.seed(task$algorithm_seed)
  # The 3:1 A/B cells use the ordinary-Louvain graph reconstruction from the main refits.
  # The two stress grids and 5:1 A/B cells use the generator's igraph object.
  graph <- if (task$imbalance_ratio == 3L && task$mode %in% c("A", "B")) {
    graph_from_adjacency_matrix(d$A, mode = "undirected", diag = FALSE)
  } else {
    d$graph
  }
  louvain <- membership(cluster_louvain(graph, resolution = 1))
  set.seed(task$algorithm_seed)
  known_k <- stats::kmeans(d$X, centers = task$K_true, nstart = 10L, iter.max = 10L)
  rbind(
    score_fit(task, "Proposed", d$z, full$partition, full$iterations_run, full$stop_reason),
    score_fit(task, "Pure_aSBM", d$z, pure$partition, pure$iterations, pure$stop_reason),
    score_fit(task, "Louvain", d$z, louvain),
    score_fit(task, "KMeans_knownK", d$z, known_k$cluster, known_k$iter)
  )
}

rows <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i, , drop = FALSE]
  cat(sprintf("[%d/%d] %s\n", i, nrow(tasks), task$case_id))
  metrics <- run_case(task)
  write.csv(metrics, file.path(out, "runs", paste0(task$case_id, ".csv")), row.names = FALSE)
  metrics
})
raw <- do.call(rbind, rows)
write.csv(raw, file.path(out, "raw_metrics.csv"), row.names = FALSE)
summary <- aggregate(
  ARI ~ imbalance_ratio + mode + scenario + r + sigma + method, raw,
  function(x) c(n = length(x), ARI_mean = mean(x), ARI_sd = sd(x))
)
write.csv(cbind(summary[setdiff(names(summary), "ARI")], summary$ARI),
          file.path(out, "condition_summary.csv"), row.names = FALSE)
