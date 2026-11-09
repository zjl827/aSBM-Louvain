#!/usr/bin/env Rscript
# Reviewer 1, Major 1: finite scalar-potential diagnostics on balanced inputs.
suppressPackageStartupMessages({
  library(igraph)
  library(Matrix)
  library(mclust)
  library(mvtnorm)
})

script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
suite <- dirname(script)
code <- normalizePath(file.path(suite, "../../code"))
model_source <- file.path(dirname(dirname(code)), "asbm_louvain.R")
source(model_source, keep.source = TRUE)
source(file.path(code, "pure_asbm_vectorized.R"))
source(file.path(code, "simulation_design.R"))

# Instrument the fitting routine without redefining its shared helpers.
replace_one <- function(text, needle, replacement) {
  hits <- gregexpr(needle, text, fixed = TRUE)[[1L]]
  if (length(hits) != 1L || hits[1L] < 0L) stop("Unexpected wrapper structure: ", needle)
  sub(needle, replacement, text, fixed = TRUE)
}
fit_source <- paste(as.character(getSrcref(aSBM_Louvain_revised)), collapse = "\n")
observer <- paste0("aSBM_Louvain_traced <- ", fit_source)
observer <- replace_one(observer, "    verbose = FALSE) {", "    verbose = FALSE, diagnostic_alpha = 1) {")
observer <- replace_one(observer, "  for (t in seq_len(T_outer)) {",
                        "  diagnostic_previous <- NULL\n  diagnostic_rows <- list()\n  for (t in seq_len(T_outer)) {")
needle <- "    if (length(elbo_history) == 0L || K_new != tail(K_history, 1L)) {"
fragment <- paste(readLines(file.path(suite, "observer_fragment.R")), collapse = "\n")
observer <- replace_one(observer, needle, paste0(fragment, "\n\n", needle))
observer <- replace_one(observer, "    K_history = K_history\n",
                        "    K_history = K_history,\n    diagnostics = do.call(rbind, diagnostic_rows)\n")
eval(parse(text = observer))

args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript experiments/potential_diagnostics/run.R OUTPUT_DIR [CASE_ID ...]")
out <- args[1L]
tasks <- read.csv(file.path(suite, "../main_simulations/tasks.csv"))
stopifnot(nrow(tasks) == 480L, !anyDuplicated(tasks$case_id))
if (length(args) > 1L) {
  stopifnot(all(args[-1L] %in% tasks$case_id))
  tasks <- tasks[tasks$case_id %in% args[-1L], ]
}
dir.create(file.path(out, "runs"), recursive = TRUE, showWarnings = FALSE)

run_case <- function(task) {
  d <- simulate_asbm_data_revised(task$N, task$K_true, 50L, task$r, task$sigma, 20, task$data_seed)
  set.seed(task$algorithm_seed)
  initial <- pure_asbm_matched_kmeans_initialization(d$A, d$X, resolution = 1, nstart = 5L)$partition
  post_init_rng <- .Random.seed
  options <- list(d$A, d$X, T_outer = 1000L, rho_tau = 0.8, tau0 = 1, resolution = 1,
                  epsilon_Q = 1e-4, epsilon_L = 1e-3, k_stable_required = 3L,
                  init_method = "kmeans", initial_partition = initial, prior_feedback = TRUE)
  original <- do.call(aSBM_Louvain_revised, options)
  final_rng <- .Random.seed
  assign(".Random.seed", post_init_rng, envir = .GlobalEnv)
  traced <- do.call(aSBM_Louvain_traced, c(options, list(diagnostic_alpha = 1)))
  stopifnot(identical(final_rng, .Random.seed),
            isTRUE(all.equal(original, traced[names(original)], tolerance = 1e-12)))
  diagnostic <- traced$diagnostics
  diagnostic$case_id <- task$case_id
  stopifnot(all(diagnostic$relative_identity_residual[diagnostic$eligible] < 1e-8))
  diagnostic
}

rows <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i, , drop = FALSE]
  cat(sprintf("[%d/%d] %s\n", i, nrow(tasks), task$case_id))
  diagnostic <- run_case(task)
  write.csv(diagnostic, file.path(out, "runs", paste0(task$case_id, ".csv")), row.names = FALSE)
  diagnostic
})
trace <- do.call(rbind, rows)
write.csv(trace, file.path(out, "iteration_diagnostics.csv"), row.names = FALSE)
eligible <- trace[trace$eligible, ]
components <- c("I", "II", "III", "IV", "delta_potential")
counts <- lapply(components, function(component) {
  decreased <- eligible[[component]] < -1e-8
  data.frame(component = component, cases = nrow(tasks), eligible_transitions = nrow(eligible),
             decrease_count = sum(decreased), runs_with_decrease = length(unique(eligible$case_id[decreased])))
})
write.csv(do.call(rbind, counts), file.path(out, "potential_summary.csv"), row.names = FALSE)
