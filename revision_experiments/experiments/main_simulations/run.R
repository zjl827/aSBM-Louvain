#!/usr/bin/env Rscript
# Refit the four methods on each independently generated balanced network.
suppressPackageStartupMessages({
  source("code/asbmlou.R")
  source("code/asbmlou_revised_stop.R")
  source("code/simulation_design.R")
})

source("code/pure_asbm_vectorized.R")

tasks <- read.csv("experiments/main_simulations/tasks.csv", stringsAsFactors = FALSE)
case_ids <- commandArgs(trailingOnly = TRUE)
if (length(case_ids)) {
  if (!all(case_ids %in% tasks$case_id)) stop("Unknown case ID in command line.")
  tasks <- tasks[tasks$case_id %in% case_ids, , drop = FALSE]
}

out <- Sys.getenv("MAIN_SIM_OUTPUT", "outputs/main_simulations")
dir.create(file.path(out, "runs"), recursive = TRUE, showWarnings = FALSE)

run_case <- function(task) {
  path <- file.path(out, "runs", paste0(task$case_id, ".rds"))
  if (file.exists(path)) {
    saved <- readRDS(path)
    stopifnot(identical(saved$task, task))
    return(saved$metrics)
  }

  dat <- simulate_asbm_data_revised(
    task$N, task$K_true, 50L, task$r, task$sigma, 20, task$data_seed
  )
  seed <- task$algorithm_seed
  set.seed(seed)
  initial <- pure_asbm_matched_kmeans_initialization(dat$A, dat$X)$partition
  post_init_rng <- .Random.seed
  set.seed(seed)
  graph <- igraph::graph_from_adjacency_matrix(dat$A, mode = "undirected", diag = FALSE)
  louvain <- igraph::membership(igraph::cluster_louvain(graph, resolution = 1))
  set.seed(seed)
  kmeans <- stats::kmeans(dat$X, centers = task$K_true, nstart = 10L, iter.max = 10L)

  assign(".Random.seed", post_init_rng, envir = .GlobalEnv)
  proposed <- aSBM_Louvain_revised(dat$A, dat$X, initial_partition = initial)
  assign(".Random.seed", post_init_rng, envir = .GlobalEnv)
  pure <- fit_pure_asbm_vectorized(
    dat$A, dat$X, initial_partition = initial, tolerance = 1,
    timeout_sec = 3600, max_iterations = 1000L
  )

  partitions <- list(Proposed = proposed$partition, KMeans_knownK = kmeans$cluster,
                     Louvain = louvain, Pure_aSBM = pure$partition)
  metrics <- cbind(
    task[rep(1L, length(partitions)), ],
    method = names(partitions),
    ARI = vapply(partitions, function(z) mclust::adjustedRandIndex(dat$z, z), numeric(1)),
    final_K = lengths(lapply(partitions, unique)),
    iterations = c(proposed$iterations_run, kmeans$iter, NA_integer_, pure$iterations),
    stop_reason = c(proposed$stop_reason, "reference_partition", "reference_partition",
                    pure$stop_reason)
  )
  tmp <- tempfile(pattern = paste0(task$case_id, "-"), tmpdir = dirname(path))
  saveRDS(list(task = task, metrics = metrics, partitions = partitions,
               initial_partition = initial, post_init_rng = post_init_rng), tmp)
  stopifnot(file.rename(tmp, path))
  message("Completed ", task$case_id)
  metrics
}

rows <- lapply(seq_len(nrow(tasks)), function(i) run_case(tasks[i, , drop = FALSE]))
filename <- if (length(case_ids)) "selected_metrics.csv" else "raw_metrics.csv"
write.csv(do.call(rbind, rows), file.path(out, filename), row.names = FALSE)
