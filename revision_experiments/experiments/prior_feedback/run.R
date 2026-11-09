#!/usr/bin/env Rscript
# Supplement S7: paired Full versus Uniform prior, with matched inputs and RNG.
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
source(file.path(code, "asbmlou_revised_stop.R"))
source(file.path(code, "pure_asbm_vectorized.R"))
source(file.path(code, "simulation_design.R"))

args <- commandArgs(TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript experiments/prior_feedback/run.R OUTPUT_DIR {balanced|real|all} [CASE_ID ...]")
}
out <- args[1L]
scope <- match.arg(args[2L], c("balanced", "real", "all"))
requested <- args[-c(1L, 2L)]
real_dir <- Sys.getenv("REVISION_REALDATA_DIR", "results/realdata_current/inputs")
orientation <- "expression_graph_spatial_attributes"

balanced <- read.csv(file.path(suite, "../main_simulations/tasks.csv"))
stopifnot(nrow(balanced) == 480L, !anyDuplicated(balanced$case_id))
tasks <- transform(balanced[, c("case_id", "scenario", "point_id", "N", "K_true", "r",
                                  "sigma", "rep", "data_seed", "algorithm_seed")],
                   scope = "balanced", dataset = "synthetic", resolution = 1)
for (dataset in c("STARmap", "MERFISH")) {
  runs <- seq_len(10L)
  tasks <- rbind(tasks, data.frame(
    case_id = paste(dataset, orientation, sprintf("rep%02d", runs), sep = "__"),
    scenario = dataset, point_id = dataset, N = if (dataset == "STARmap") 1207L else 5803L,
    K_true = NA_integer_, r = NA_real_, sigma = NA_real_, rep = runs,
    data_seed = NA_integer_, algorithm_seed = 20260901L + runs, scope = "real",
    dataset = dataset, resolution = if (dataset == "STARmap") 1.5 else 0.5
  ))
}
if (scope != "all") tasks <- tasks[tasks$scope == scope, ]
if (length(requested)) {
  stopifnot(all(requested %in% tasks$case_id))
  tasks <- tasks[tasks$case_id %in% requested, ]
}
inputs <- lapply(unique(tasks$dataset[tasks$scope == "real"]), function(dataset) {
  path <- file.path(real_dir, paste0(dataset, "__", orientation, ".rds"))
  if (!file.exists(path)) stop("Missing input: ", path, ". Run experiments/realdata_current/prepare.R.")
  readRDS(path)
})
names(inputs) <- unique(tasks$dataset[tasks$scope == "real"])
dir.create(file.path(out, "runs"), recursive = TRUE, showWarnings = FALSE)

run_pair <- function(task) {
  if (task$scope == "balanced") {
    data <- simulate_asbm_data_revised(task$N, task$K_true, 50L, task$r, task$sigma, 20, task$data_seed)
    A <- data$A
    X <- data$X
    truth <- data$z
  } else {
    data <- inputs[[task$dataset]]
    A <- data$adjacency
    X <- data$attributes
    truth <- data$truth
  }
  stopifnot(nrow(A) == ncol(A), nrow(A) == nrow(X), length(truth) == task$N, nrow(X) == task$N)
  set.seed(task$algorithm_seed)
  init <- pure_asbm_matched_kmeans_initialization(A, X, resolution = task$resolution, nstart = 5L)
  post_init_rng <- .Random.seed
  fit <- function(prior_feedback) {
    assign(".Random.seed", post_init_rng, envir = .GlobalEnv)
    aSBM_Louvain_revised(
      A, X, T_outer = 1000L, rho_tau = 0.8, tau0 = 1, resolution = task$resolution,
      epsilon_Q = 1e-4, epsilon_L = 1e-3, k_stable_required = 3L,
      init_method = "kmeans", initial_partition = init$partition, prior_feedback = prior_feedback
    )
  }
  full <- fit(TRUE)
  uniform <- fit(FALSE)
  stopifnot(identical(full$initial_partition, uniform$initial_partition))
  score <- function(fit) {
    c(ARI = adjustedRandIndex(truth, fit$partition),
      NMI = compare(as.integer(factor(truth)), as.integer(factor(fit$partition)), method = "nmi"))
  }
  a <- score(full)
  b <- score(uniform)
  data.frame(task, initial_K = init$estimated_groups, initial_ARI = adjustedRandIndex(truth, init$partition),
             full_ARI = a[["ARI"]], uniform_ARI = b[["ARI"]], ARI_delta_full_minus_uniform = a[["ARI"]] - b[["ARI"]],
             full_NMI = a[["NMI"]], uniform_NMI = b[["NMI"]], NMI_delta_full_minus_uniform = a[["NMI"]] - b[["NMI"]],
             full_K = full$final_K, uniform_K = uniform$final_K,
             full_iterations = full$iterations_run, uniform_iterations = uniform$iterations_run,
             full_stop_reason = full$stop_reason, uniform_stop_reason = uniform$stop_reason)
}

rows <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i, , drop = FALSE]
  cat(sprintf("[%d/%d] %s\n", i, nrow(tasks), task$case_id))
  metrics <- run_pair(task)
  write.csv(metrics, file.path(out, "runs", paste0(task$case_id, ".csv")), row.names = FALSE)
  metrics
})
raw <- do.call(rbind, rows)
write.csv(raw, file.path(out, "paired_raw.csv"), row.names = FALSE)
groups <- split(raw, raw$point_id)
summary <- lapply(groups, function(x) {
  columns <- c("full_ARI", "uniform_ARI", "ARI_delta_full_minus_uniform", "full_NMI", "uniform_NMI",
               "NMI_delta_full_minus_uniform")
  averages <- setNames(lapply(x[columns], mean), paste0(columns, "_mean"))
  deviations <- setNames(lapply(x[columns], sd), paste0(columns, "_sd"))
  data.frame(x[1L, c("scope", "dataset", "scenario", "point_id", "r", "sigma")],
             n = nrow(x), averages, deviations)
})
write.csv(do.call(rbind, summary), file.path(out, "pointwise_summary.csv"), row.names = FALSE)
