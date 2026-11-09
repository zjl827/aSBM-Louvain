#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(igraph); library(mvtnorm); library(mclust)})
args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[startsWith(args, "--file=")][1]))
script_dir <- dirname(script)
workspace <- normalizePath(file.path(script_dir, "..", ".."))
code <- file.path(workspace, "code")
source(file.path(code, "asbmlou.R"))
source(file.path(code, "asbmlou_revised_stop.R"))
native <- file.path(code, paste0("pure_asbm_floor_correction", .Platform$dynlib.ext))
if (!file.exists(native)) stop("Compile the aSBM baseline with R CMD SHLIB code/pure_asbm_floor_correction.c")
if (!is.loaded("pure_asbm_floor_correction_c")) dyn.load(native)
source(file.path(code, "pure_asbm_vectorized.R"))
source(file.path(workspace, "experiments", "weight_misspecification", "design_helpers.R"))
stopifnot(is.loaded("pure_asbm_floor_correction_c"))

out <- Sys.getenv("PILOT_OUTPUT", file.path(workspace, "outputs", "joint_signals"))
dir.create(file.path(out, "runs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "data_cache"), recursive = TRUE, showWarnings = FALSE)
grid <- read.csv(file.path(script_dir, "grid.csv"),
                 stringsAsFactors = FALSE)
reps <- as.integer(Sys.getenv("PILOT_REPS", "10"))
workers <- as.integer(Sys.getenv("PILOT_WORKERS", "4"))
smoke <- Sys.getenv("PILOT_SMOKE", "0") == "1"
if (smoke) { grid <- grid[grid$condition_id %in% c("C01", "C14", "C16"), ]; reps <- 1L }
r_values <- c(1.5, 3, 4, 7.5)
version <- "proposed_consistent_pilot_v1"
settings <- list(stage = if (smoke) "smoke" else "pilot", repetitions = reps,
                 workers = workers, version = version, N = 1000L, K = 5L, D = 50L,
                 tau0 = 1, rho_tau = 0.8, resolution = 1,
                 epsilon_Q = 1e-4, epsilon_L = 1e-3,
                 k_stable_required = 3L, outer_cap = 1000L,
                 effective_tau_floor = 1e-5, pure_tolerance = 1,
                 per_fit_timeout_sec = 300, RNG_kind = RNGkind(),
                 R_session = capture.output(sessionInfo()),
                 runtime_warning = "Parallel diagnostic timings; native-accelerated Pure baseline; not benchmark timings.")
saveRDS(settings, file.path(out, "run_settings.rds"))
capture.output(str(settings), file = file.path(out, "run_settings.txt"))

fit_safely <- function(fun) {
  warning_messages <- character()
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch({
    setTimeLimit(cpu = Inf, elapsed = 300, transient = TRUE)
    withCallingHandlers(fun(), warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  }, error = function(e) structure(list(message = conditionMessage(e)), class = "pilot_error"),
  finally = setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  list(fit = fit, elapsed = proc.time()[["elapsed"]] - start,
       warnings = paste(unique(warning_messages), collapse = " | "))
}

run_group <- function(task) {
  ri <- task$ri; rep_id <- task$rep; ratio <- r_values[ri]
  conditions <- grid[grid$r == ratio, ]
  n <- 1000L; k <- 5L; d <- 50L; sizes <- rep(200L, k)
  truth <- rep(seq_len(k), times = sizes)
  seeds <- list(graph_seed = 92100000L + ri * 100000L + rep_id,
                means_seed = 92101000L + rep_id,
                residual_seed = 92102000L + rep_id,
                assignment_seed = 92103000L + rep_id,
                graph_algorithm_seed = 92500000L + ri * 1000L + rep_id,
                algorithm_seed = 92600000L + rep_id)
  set.seed(seeds$graph_seed)
  p_out <- 20 / (199 * ratio + 800)
  p_in <- ratio * p_out
  block <- matrix(p_out, k, k); diag(block) <- p_in
  g <- igraph::sample_sbm(n, pref.matrix = block, block.sizes = sizes,
                          directed = FALSE, loops = FALSE)
  A <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  set.seed(seeds$means_seed)
  means <- matrix(rnorm(k * d), k, d)
  set.seed(seeds$residual_seed)
  residual <- matrix(rnorm(n * d), n, d)
  set.seed(seeds$assignment_seed)
  high <- sort(sample(seq_len(k), 2L, replace = FALSE))
  set.seed(seeds$graph_algorithm_seed)
  start <- proc.time()[["elapsed"]]
  lp <- as.integer(factor(igraph::membership(igraph::cluster_louvain(g, resolution = 1))))
  graph_init_elapsed <- proc.time()[["elapsed"]] - start
  initial_K <- max(2L, length(unique(lp)))
  graph_hash <- hash_object(A)
  saveRDS(list(A = A, means = means, residual = residual, high_groups = high,
               truth = truth, seeds = seeds, r = ratio, p_in = p_in, p_out = p_out),
          file.path(out, "data_cache", sprintf("r%d_rep%02d.rds", ri, rep_id)))

  completed <- character()
  for (row_index in seq_len(nrow(conditions))) {
    cnd <- conditions[row_index, ]
    key <- sprintf("%s_rep%02d", cnd$condition_id, rep_id)
    destination <- file.path(out, "runs", paste0(key, ".rds"))
    if (file.exists(destination)) {
      old <- readRDS(destination)
      stopifnot(old$version == version, nrow(old$metrics) == 4L)
      completed <- c(completed, key)
      next
    }
    sigmas <- community_noise_sd(cnd$high_low_sd_ratio, cnd$sigma_rms, high, sizes)
    X <- means[truth, , drop = FALSE] + residual * sigmas[truth]
    stopifnot(abs(mean(sigmas^2) - cnd$sigma_rms^2) < 1e-10,
              max(abs(A - t(A))) == 0, all(diag(A) == 0))
    set.seed(seeds$algorithm_seed)
    start <- proc.time()[["elapsed"]]
    km_attempt <- fit_safely(function() stats::kmeans(X, centers = initial_K, nstart = 5))
    init_elapsed <- graph_init_elapsed + proc.time()[["elapsed"]] - start
    km_failed <- inherits(km_attempt$fit, "pilot_error")
    initial <- if (km_failed) NULL else as.integer(factor(km_attempt$fit$cluster))
    post_rng <- .Random.seed
    initial_hash <- if (km_failed) NA_character_ else hash_object(initial)
    attribute_hash <- hash_object(X)
    meta <- cnd
    meta$rep <- rep_id
    meta$initial_K <- initial_K
    meta$high_groups <- paste(high, collapse = ":")
    meta$noise_sd <- paste(sprintf("%.12g", sigmas), collapse = ":")
    meta$graph_hash <- graph_hash; meta$attribute_hash <- attribute_hash
    meta$initial_partition_hash <- initial_hash
    meta$algorithm_version <- version
    meta$initialization_elapsed_sec <- init_elapsed
    meta$realized_mean_degree <- sum(A) / n
    meta$p_in <- p_in; meta$p_out <- p_out
    for (nm in names(seeds)) meta[[nm]] <- seeds[[nm]]
    meta$post_init_rng_hash <- hash_object(post_rng)
    outputs <- list(); local <- list(); predictions <- list()
    methods <- c("Proposed", "Pure_aSBM", "Louvain", "KMeans_estK")
    for (method in methods) {
      if (method == "Louvain") {
        attempt <- list(fit = list(partition = lp, final_K = max(lp), iterations = NA_integer_,
                                   stop_reason = "reference_partition"),
                        elapsed = graph_init_elapsed, warnings = "")
      } else if (km_failed) {
        attempt <- list(fit = km_attempt$fit, elapsed = NA_real_, warnings = km_attempt$warnings)
      } else if (method == "KMeans_estK") {
        attempt <- list(fit = list(partition = initial, final_K = initial_K,
                                   iterations = km_attempt$fit$iter, stop_reason = "reference_partition"),
                        elapsed = km_attempt$elapsed, warnings = km_attempt$warnings)
      } else {
        assign(".Random.seed", post_rng, envir = .GlobalEnv)
        if (method == "Proposed") {
          attempt <- fit_safely(function() aSBM_Louvain_revised(
            A, X, T_outer = 1000L, tau0 = 1, rho_tau = 0.8, resolution = 1,
            epsilon_Q = 1e-4, epsilon_L = 1e-3, k_stable_required = 3L,
            init_method = "kmeans", initial_partition = initial, verbose = FALSE))
        } else {
          attempt <- fit_safely(function() fit_pure_asbm_vectorized(
            A, X, initial_partition = initial, tolerance = 1,
            timeout_sec = 300, max_iterations = 1000L))
        }
      }
      failed <- inherits(attempt$fit, "pilot_error")
      metrics <- meta
      metrics$method <- method
      metrics$elapsed_sec <- attempt$elapsed
      metrics$warnings <- attempt$warnings
      metrics$error <- if (failed) attempt$fit$message else ""
      metrics$status <- if (failed) "error" else "ok"
      if (failed) {
        metrics$ARI <- metrics$NMI <- NA_real_
        metrics$final_K <- metrics$iterations <- NA_integer_
        metrics$stop_reason <- "error"
        metrics$met_stopping_rule <- FALSE
        comm <- data.frame(community_id = seq_len(k), precision = NA_real_, recall = NA_real_,
                           F1 = NA_real_, is_high_noise = seq_len(k) %in% high,
                           true_size = sizes, predicted_size = NA_integer_, matched_label = NA_integer_)
      } else {
        fit <- attempt$fit
        pred <- as.integer(fit$partition)
        stopifnot(length(pred) == n, !anyNA(pred))
        predictions[[method]] <- pred
        metrics$ARI <- mclust::adjustedRandIndex(truth, pred)
        metrics$NMI <- partition_nmi(truth, pred)
        metrics$final_K <- length(unique(pred))
        metrics$iterations <- if (method == "Proposed") fit$iterations_run else fit$iterations
        metrics$stop_reason <- fit$stop_reason
        metrics$met_stopping_rule <- if (method == "Proposed") fit$stop_reason == "joint_stability" else if (method == "Pure_aSBM") fit$stop_reason %in% c("tolerance", "objective_not_increased") else NA
        comm <- community_scores(truth, pred, high)
      }
      comm$condition_id <- cnd$condition_id; comm$rep <- rep_id; comm$method <- method
      outputs[[method]] <- metrics; local[[method]] <- comm
    }
    result <- list(version = version, metrics = do.call(rbind, outputs),
                   communities = do.call(rbind, local), predictions = predictions,
                   initial_partition = initial, post_init_rng = post_rng,
                   high_groups = high, sigmas = sigmas, seeds = seeds)
    tmp <- paste0(destination, ".tmp")
    saveRDS(result, tmp)
    stopifnot(file.rename(tmp, destination))
    write.csv(result$metrics, file.path(out, "runs", paste0(key, "_metrics.csv")), row.names = FALSE)
    cat(sprintf("DONE %s r=%g sigma=%g h=%g Proposed=%.3f Pure=%.3f errors=%d\n",
                key, ratio, cnd$sigma_rms, cnd$high_low_sd_ratio,
                result$metrics$ARI[result$metrics$method == "Proposed"],
                result$metrics$ARI[result$metrics$method == "Pure_aSBM"],
                sum(result$metrics$status != "ok")))
    flush.console()
    completed <- c(completed, key)
  }
  completed
}

tasks <- list()
for (ri in seq_along(r_values)) if (r_values[ri] %in% grid$r) {
  for (rep_id in seq_len(reps)) tasks[[length(tasks) + 1L]] <- list(ri = ri, rep = rep_id)
}
cat(sprintf("START %s: conditions=%d repeats=%d workers=%d\n", version, nrow(grid), reps, workers))
flush.console()
results <- parallel::mclapply(tasks, run_group, mc.cores = workers,
                             mc.preschedule = FALSE, mc.set.seed = FALSE)
if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
  stop("Worker failure; completed cases retained for restart.")
}
paths <- list.files(file.path(out, "runs"), pattern = "\\.rds$", full.names = TRUE)
all_runs <- lapply(paths, readRDS)
metrics <- do.call(rbind, lapply(all_runs, `[[`, "metrics"))
communities <- do.call(rbind, lapply(all_runs, `[[`, "communities"))
stopifnot(nrow(metrics) == nrow(grid) * reps * 4L,
          !anyDuplicated(paste(metrics$condition_id, metrics$rep, metrics$method)))
write.csv(metrics, file.path(out, "raw_metrics.csv"), row.names = FALSE)
write.csv(communities, file.path(out, "community_metrics.csv"), row.names = FALSE)
saveRDS(lapply(all_runs, function(x) x[c("predictions", "initial_partition", "high_groups", "sigmas", "seeds")]),
        file.path(out, "predictions_and_initializations.rds"))
cat(sprintf("COMPLETE rows=%d communities=%d errors=%d\n", nrow(metrics), nrow(communities), sum(metrics$status != "ok")))
