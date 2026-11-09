#!/usr/bin/env Rscript
# Refit the four graph/attribute methods used in the orientation experiment.
suppressPackageStartupMessages({library(igraph); library(Matrix); library(mclust); library(mvtnorm)})
algorithm <- "../asbm_louvain.R"
source(algorithm)
source("code/pure_asbm_vectorized.R")
stopifnot(formals(stats::kmeans)$iter.max == 10L, eval(formals(stats::kmeans)$algorithm)[1L] == "Hartigan-Wong")
input <- Sys.getenv("REVISION_REALDATA_DIR", "results/realdata_current/inputs")
output <- Sys.getenv("REVISION_REALDATA_OUT", "results/realdata_current")
args <- commandArgs(TRUE)
stopifnot(!length(args) || identical(args, "smoke"))
runs <- if (length(args)) 1L else seq_len(10L)
datasets <- if (length(args)) "STARmap" else c("STARmap", "MERFISH")
orientations <- c("expression_graph_spatial_attributes", "spatial_graph_expression_attributes")
if (length(args)) orientations <- orientations[1L]
methods <- c("Proposed", "Pure_aSBM", "Louvain", "KMeans_estK")
dir.create(file.path(output, "runs"), showWarnings = FALSE, recursive = TRUE)
core_files <- c(algorithm, "code/pure_asbm_vectorized.R",
  "experiments/realdata_current/run.R")
code_md5 <- tools::md5sum(core_files)
capture.output(sessionInfo(), file = file.path(output, "R_sessionInfo.txt"))

capture_fit <- function(fn) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(fn(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) structure(list(message = conditionMessage(e)), class = "fit_error"))
  list(value = value, warnings = paste(unique(warnings), collapse = " | "))
}

run_case <- function(dat, dataset, orientation, rep) {
  A <- dat$adjacency; X <- dat$attributes; truth <- dat$truth
  seed <- 20260901L + rep
  set.seed(seed)
  init <- capture_fit(function() pure_asbm_matched_kmeans_initialization(A, X, dat$resolution, nstart = 5L))
  rng <- .Random.seed
  initial <- init$value$partition
  # Keep the independent originally submitted Louvain protocol, including cell order.
  ord <- order(dat$cell_ids, method = "radix")
  fits <- lapply(methods, function(method) {
    if (method == "Louvain") {
      return(capture_fit(function() {
        g <- graph_from_adjacency_matrix(A[ord, ord], mode = "undirected", diag = FALSE)
        set.seed(123L + rep)
        z <- integer(length(ord))
        z[ord] <- as.integer(membership(cluster_louvain(g, resolution = 1)))
        list(partition = z, stop_reason = "reference_partition")
      }))
    }
    if (is.null(initial)) return(init)
    assign(".Random.seed", rng, envir = .GlobalEnv)
    if (method == "KMeans_estK") {
      return(list(value = list(partition = initial, stop_reason = "reference_partition"), warnings = ""))
    }
    capture_fit(function() {
      if (method == "Proposed") {
        aSBM_Louvain_revised(A, X, T_outer = 1000L, rho_tau = .8, tau0 = 1, resolution = dat$resolution,
          epsilon_Q = 1e-4, epsilon_L = 1e-3, k_stable_required = 3L,
          init_method = "kmeans", initial_partition = initial, prior_feedback = TRUE, verbose = FALSE)
      } else {
        fit_pure_asbm_vectorized(A, X, initial_partition = initial, tolerance = 1,
          timeout_sec = 3600, max_iterations = 1000L)
      }
    })
  })
  names(fits) <- methods
  metrics <- do.call(rbind, lapply(methods, function(method) {
    f <- fits[[method]]$value
    fail <- inherits(f, "fit_error")
    z <- f$partition
    stop_reason <- if (fail) "error" else f$stop_reason
    iterations <- if (method == "Proposed") f$iterations_run else f$iterations
    data.frame(dataset = dataset, orientation = orientation, rep = rep, method = method,
      seed = if (method == "Louvain") 123L + rep else seed,
      resolution = if (method == "Louvain") 1 else dat$resolution,
      N = nrow(X), D = ncol(X),
      ARI = if (fail) NA_real_ else adjustedRandIndex(truth, z),
      NMI = if (fail) NA_real_ else igraph::compare(as.numeric(factor(truth)), as.numeric(factor(z)), method = "nmi"),
      initial_K = if (method == "Louvain" || is.null(initial)) NA_integer_ else length(unique(initial)),
      final_K = if (fail) NA_integer_ else length(unique(z)),
      iterations = if (is.null(iterations)) NA_integer_ else iterations, stop_reason = stop_reason,
      status = if (fail) "error" else if (stop_reason %in% c("safety_cap", "iteration_cap", "timeout")) "incomplete" else "ok",
      init_warnings = if (method == "Louvain") "" else init$warnings, warnings = fits[[method]]$warnings,
      error = if (fail) f$message else "")
  }))
  list(metrics = metrics, initial_partition = initial, post_init_rng = rng, fits = fits, cell_ids = dat$cell_ids)
}

rows <- list()
for (dataset in datasets) for (orientation in orientations) {
  key <- paste(dataset, orientation, sep = "__")
  path <- file.path(input, paste0(key, ".rds"))
  dat <- readRDS(path)
  stopifnot(length(dat$cell_ids) == nrow(dat$attributes), !anyDuplicated(dat$cell_ids))
  input_md5 <- unname(tools::md5sum(path))
  for (rep in runs) {
    id <- paste(key, sprintf("rep%02d", rep), sep = "__")
    dest <- file.path(output, "runs", paste0(id, ".rds"))
    if (file.exists(dest)) {
      record <- readRDS(dest)
      if (!identical(record$input_md5, input_md5) || !identical(record$code_md5, code_md5)) {
        stop("Saved real-data checkpoints use different code or inputs. Set REVISION_REALDATA_OUT to a fresh directory.")
      }
    } else {
      record <- run_case(dat, dataset, orientation, rep)
      record$input_md5 <- input_md5
      record$code_md5 <- code_md5
      saveRDS(record, paste0(dest, ".tmp"))
      stopifnot(file.rename(paste0(dest, ".tmp"), dest))
    }
    rows[[id]] <- record$metrics
    message("Finished ", id)
  }
}
raw <- do.call(rbind, rows)
write.csv(raw, file.path(output, "raw_metrics.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(raw, interaction(raw$dataset, raw$orientation, raw$method, drop = TRUE)), function(d) {
  data.frame(d[1L, c("dataset", "orientation", "method")],
    n = nrow(d), n_returned = sum(is.finite(d$ARI)), n_completed = sum(d$status == "ok"),
    ARI_mean = mean(d$ARI, na.rm = TRUE), ARI_SD = sd(d$ARI, na.rm = TRUE),
    NMI_mean = mean(d$NMI, na.rm = TRUE), NMI_SD = sd(d$NMI, na.rm = TRUE))
}))
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
