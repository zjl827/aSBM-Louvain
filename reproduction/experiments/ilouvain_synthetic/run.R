#!/usr/bin/env Rscript
# Reviewer 1 and Supplement S9: I-Louvain on the original synthetic generator.
suppressPackageStartupMessages({
  library(igraph)
  library(mclust)
})

script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]))
suite <- dirname(script)
source(file.path(suite, "../../code/simulation_design.R"))

args <- commandArgs(TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript experiments/ilouvain_synthetic/run.R EXTERNAL_IL_PY OUTPUT_DIR [CASE_ID ...]")
}
il_path <- normalizePath(args[1L], mustWork = TRUE)
out <- args[2L]
python <- Sys.getenv("PYTHON", "python3")

columns <- c("case_id", "scenario", "N", "K_true", "r", "sigma", "rep", "data_seed")
balanced <- read.csv(file.path(suite, "../main_simulations/tasks.csv"))
unequal <- read.csv(file.path(suite, "../unequal_sizes/tasks.csv"))
tasks <- rbind(
  data.frame(balanced[columns], scope = "balanced", block_sizes = NA_character_, array_attributes = FALSE),
  data.frame(unequal[columns], scope = "unequal", block_sizes = unequal$block_sizes,
             array_attributes = unequal$imbalance_ratio == 5L | unequal$mode == "C")
)
stopifnot(nrow(tasks) == 840L, !anyDuplicated(tasks$case_id))
requested <- args[-c(1L, 2L)]
if (length(requested)) {
  stopifnot(all(requested %in% tasks$case_id))
  tasks <- tasks[tasks$case_id %in% requested, ]
}
dir.create(file.path(out, "runs"), recursive = TRUE, showWarnings = FALSE)

run_case <- function(task) {
  if (task$scope == "balanced") {
    d <- simulate_asbm_data_revised(task$N, task$K_true, 50L, task$r, task$sigma, 20, task$data_seed)
  } else {
    sizes <- as.integer(strsplit(task$block_sizes, ":", fixed = TRUE)[[1L]])
    d <- simulate_unbalanced_data(sizes, 50L, task$r, task$sigma, 20, task$data_seed)
  }
  temp <- tempfile("ilouvain_case_")
  dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE))
  edges <- file.path(temp, "edges.bin")
  attributes <- file.path(temp, "attributes.bin")
  labels <- file.path(temp, "labels.bin")
  writeBin(as.integer(t(as_edgelist(d$graph, names = FALSE) - 1L)), edges, size = 4L, endian = "little")
  writeBin(as.double(t(d$X)), attributes, size = 8L, endian = "little")

  command <- c(file.path(suite, "run_case.py"), "--ilouvain-code", il_path,
               "--edges", edges, "--attributes", attributes, "--labels-out", labels,
               "--n", task$N, "--d", ncol(d$X))
  if (task$array_attributes) command <- c(command, "--array-attributes")
  log <- system2(python, shQuote(command), stdout = TRUE, stderr = TRUE, timeout = 3600)
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) {
    stop("I-Louvain failed for ", task$case_id, "\n", paste(log, collapse = "\n"))
  }
  pred <- readBin(labels, integer(), n = task$N, size = 4L, endian = "little")
  stopifnot(length(pred) == task$N, !anyNA(pred))
  info <- grep("^external_code_sha256=", log, value = TRUE)
  stopifnot(length(info) == 1L)
  data.frame(task, edges = ecount(d$graph), K_found = length(unique(pred)),
             ARI = adjustedRandIndex(d$z, pred),
             NMI = compare(as.integer(factor(d$z)), as.integer(factor(pred)), method = "nmi"),
             external_implementation = info)
}

rows <- lapply(seq_len(nrow(tasks)), function(i) {
  task <- tasks[i, , drop = FALSE]
  result <- run_case(task)
  write.csv(result, file.path(out, "runs", paste0(task$case_id, ".csv")), row.names = FALSE)
  cat(sprintf("[%d/%d] %s K=%d ARI=%.6f\n", i, nrow(tasks), task$case_id, result$K_found, result$ARI))
  result
})
write.csv(do.call(rbind, rows), file.path(out, "raw_metrics.csv"), row.names = FALSE)
