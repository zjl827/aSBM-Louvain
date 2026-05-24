#!/usr/bin/env Rscript
# ==============================================================================
# Unified simulation runner for aSBM-Louvain paper.
#
# Usage:
#   Rscript run_simulations.R baselines A   [sigma]     # A scenarios baseline
#   Rscript run_simulations.R baselines B                # B scenarios baseline
#   Rscript run_simulations.R pure      A   [sigma]     # Pure-aSBM A scenarios
#   Rscript run_simulations.R pure      B   [r]         # Pure-aSBM B scenarios
#   Rscript run_simulations.R all                       # Run everything
#
# Defaults: sigma=3.0, r=4.0, 5 reps, 4 cores, 10h timeout
# ==============================================================================

library(igraph)
if (!requireNamespace("mclust", quietly = TRUE)) stop("Package 'mclust' required.")

# ---- Data generation ----
simulate_asbm_data <- function(N, K, D, ratio, sigma, target_degree, seed) {
  set.seed(seed)
  block_sizes <- rep(N %/% K, K)
  if (N %% K > 0) block_sizes[seq_len(N %% K)] <- block_sizes[seq_len(N %% K)] + 1
  group_size <- N / K
  p_out <- target_degree / (((group_size - 1) * ratio) + (N - group_size))
  p_in <- ratio * p_out
  pref_matrix <- matrix(p_out, K, K); diag(pref_matrix) <- p_in
  graph_obj <- sample_sbm(N, pref.matrix = pref_matrix, block.sizes = block_sizes,
                          directed = FALSE, loops = FALSE)
  Adj <- as.matrix(as_adjacency_matrix(graph_obj, sparse = FALSE))
  z_true <- rep(seq_len(K), times = block_sizes)
  centers <- matrix(rnorm(K * D), nrow = K, ncol = D)
  X <- matrix(0, nrow = N, ncol = D)
  for (k in seq_len(K)) {
    ids <- which(z_true == k)
    X[ids, ] <- matrix(rnorm(length(ids) * D, mean = rep(centers[k, ], each = length(ids)),
                             sd = sigma), nrow = length(ids), byrow = FALSE)
  }
  list(A = Adj, X = X, z = z_true, graph = graph_obj)
}

safe_ari <- function(truth, pred) {
  if (is.null(pred) || length(pred) != length(truth)) return(NA_real_)
  mclust::adjustedRandIndex(truth, pred)
}

# ---- Scenario definitions ----
scenarios <- list(
  A1 = list(N = 1000L, K = 5L),
  A2 = list(N = 1000L, K = 10L),
  A3 = list(N = 5000L, K = 5L),
  B1 = list(N = 1000L, K = 5L),
  B2 = list(N = 1000L, K = 10L),
  B3 = list(N = 5000L, K = 5L)
)

A_scenarios <- c("A1", "A2", "A3")
B_scenarios <- c("B1", "B2", "B3")

r_vals_A <- c(7.5, 6.0, 4.5, 3.0, 1.5, 1.0)
r_vals_A3 <- c(7.5, 6.0, 4.5)  # A3 high-r only (low r infeasible)
sigma_vals_B <- seq(0.5, 5.0, by = 0.5)

# ---- Run baselines (Louvain / KMeans / aSBM-Louvain) ----
run_baselines <- function(mode, sigma = 3.0, reps = 5L, base_seed = 20260518) {
  source("asbmlou.R")

  if (mode == "A") {
    sc_list <- A_scenarios
    task_desc <- sprintf("sigma=%.1f", sigma)
  } else {
    sc_list <- B_scenarios
    task_desc <- "r=4.0"
  }

  all_res <- list(); idx <- 1
  for (sid in sc_list) {
    sc <- scenarios[[sid]]
    if (mode == "A") {
      r_list <- if (sid == "A3") r_vals_A3 else r_vals_A
      x_vals <- r_list; x_name <- "r"
    } else {
      x_vals <- sigma_vals_B; x_name <- "sigma"
    }

    for (xv in x_vals) {
      for (rep_id in seq_len(reps)) {
        cat(sprintf("[%s] %s=%.1f rep=%d\n", sid, x_name, xv, rep_id))
        if (mode == "A") {
          dat <- simulate_asbm_data(sc$N, sc$K, 50, xv, sigma, 20,
                                     base_seed + rep_id + as.integer(xv * 100) + sc$K * 10000)
        } else {
          dat <- simulate_asbm_data(sc$N, sc$K, 50, 4.0, xv, 20,
                                     base_seed + rep_id + as.integer(xv * 1000))
        }

        lou <- as.numeric(factor(membership(cluster_louvain(dat$graph))))
        km <- as.numeric(factor(kmeans(dat$X, centers = sc$K, nstart = 10)$cluster))
        capture.output({
          fit <- aSBM_Louvain(dat$A, dat$X, T_outer = 20,
            rho_tau = 0.8, tau0 = 1.0, resolution = 1.0, init_method = "kmeans")
        })

        all_res[[idx]] <- data.frame(
          scenario = sid, sigma = if(mode=="A") sigma else xv,
          r = if(mode=="A") xv else 4.0, rep = rep_id,
          Louvain = safe_ari(dat$z, lou),
          KMeans = safe_ari(dat$z, km),
          aSBM_Louvain = safe_ari(dat$z, fit$partition),
          stringsAsFactors = FALSE)
        idx <- idx + 1
      }
    }
  }

  res <- do.call(rbind, all_res)
  prefix <- if(mode == "A") paste0("A_baselines_sigma", sigma) else "B_baselines_r4"
  write.csv(res, paste0(prefix, "_results.csv"), row.names = FALSE)

  smry <- do.call(rbind, lapply(split(res, ~ scenario + r + sigma), function(sub) {
    data.frame(scenario = sub$scenario[1], sigma = sub$sigma[1], r = sub$r[1],
               Louvain = mean(sub$Louvain), Louvain_sd = sd(sub$Louvain),
               KMeans = mean(sub$KMeans), KMeans_sd = sd(sub$KMeans),
               aSBM_Louvain = mean(sub$aSBM_Louvain),
               aSBM_Louvain_sd = sd(sub$aSBM_Louvain), stringsAsFactors = FALSE)
  }))
  write.csv(smry, paste0(prefix, "_point_summary.csv"), row.names = FALSE)

  cat(sprintf("\n=== %s Baselines Summary ===\n", mode))
  print(smry, row.names = FALSE)
  cat(sprintf("Saved: %s_results.csv, %s_point_summary.csv\n", prefix, prefix))
}

# ---- Run Pure-aSBM ----
run_pure_asbm <- function(mode, sigma = 3.0, r_fixed = 4.0, reps = 5L,
                          base_seed = 20260518, n_cores = 4L,
                          timeout_sec = 36000) {
  library(parallel)
  source("FitAttribute.R")

  build_tasks <- function() {
    tasks <- list(); idx <- 1
    if (mode == "A") {
      for (sid in A_scenarios) {
        sc <- scenarios[[sid]]
        r_list <- if (sid == "A3") r_vals_A3 else r_vals_A
        for (rv in r_list) {
          for (rep_id in seq_len(reps)) {
            tasks[[idx]] <- list(scenario = sid, N = sc$N, K = sc$K,
                                 ratio = rv, sigma = sigma, rep = rep_id, base_seed = base_seed)
            idx <- idx + 1
          }
        }
      }
    } else {
      for (sid in B_scenarios) {
        sc <- scenarios[[sid]]
        for (sv in sigma_vals_B) {
          for (rep_id in seq_len(reps)) {
            tasks[[idx]] <- list(scenario = sid, N = sc$N, K = sc$K,
                                 ratio = r_fixed, sigma = sv, rep = rep_id, base_seed = base_seed)
            idx <- idx + 1
          }
        }
      }
    }
    tasks
  }

  worker_fn <- function(task) {
    dat <- simulate_asbm_data(task$N, task$K, 50, task$ratio, task$sigma, 20,
      task$base_seed + task$rep + as.integer(task$ratio * 100) + task$K * 10000)
    lou_init <- as.numeric(factor(membership(cluster_louvain(dat$graph))))
    K_est <- length(unique(lou_init))
    km_seed <- task$base_seed + 100000 + task$rep + as.integer(task$ratio * 100)

    # k-means init
    set.seed(km_seed)
    kmeansinit <- as.numeric(factor(kmeans(dat$X, centers = K_est, nstart = 10)$cluster))
    k <- max(kmeansinit)

    StoreMean <- list(); StoreCov <- list()
    for (kk in seq_len(k)) {
      inds <- which(kmeansinit == kk)
      if (length(inds) == 0) {
        StoreMean[[kk]] <- matrix(0, ncol = ncol(dat$X), nrow = 1)
        StoreCov[[kk]] <- matrix(0, nrow = ncol(dat$X), ncol = ncol(dat$X))
      } else {
        StoreMean[[kk]] <- matrix(colSums(matrix(dat$X[inds,], ncol=ncol(dat$X))) /
                                    length(inds), ncol = ncol(dat$X))
        StoreCov[[kk]] <- EstCovMat(matrix(dat$X[inds,], ncol = ncol(dat$X)))
      }
      if (sum(StoreCov[[kk]]) == 0) StoreCov[[kk]] <- StoreCov[[kk]] + 1e-6
    }

    ProbMat <- matrix(0.01, nrow = task$N, ncol = k)
    for (i in seq_len(length(kmeansinit))) ProbMat[i, kmeansinit[i]] <- 0.99
    Alpha <- UpdateBlockParameters(ProbMat, dat$A)
    PiMat <- apply(ProbMat, 2, function(x) sum(x)/length(x)); logPi <- log(PiMat)
    Expect <- ProbMat
    Objective <- ObjectiveFunction(dat$X, StoreMean, StoreCov, Alpha, logPi, dat$A, Expect)

    proceed <- 1; diff_val <- 2; iter <- 0; timed_out <- FALSE
    setTimeLimit(cpu = timeout_sec, elapsed = timeout_sec, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))

    fit <- tryCatch({
      while (proceed > 0 && diff_val > 1) {
        iter <- iter + 1
        for (m in seq_len(k)) {
          StoreMean[[m]] <- matrix(MuUpdate(Expect, dat$X, m))
          StoreCov[[m]] <- CovMatUpdate(StoreMean[[m]], Expect, dat$X, m)
        }
        ExpectNew <- ComputeExpectation(dat$X, StoreMean, StoreCov, Alpha, k, logPi, Expect, dat$A)
        PiMat <- colSums(ExpectNew)/nrow(ExpectNew)
        PiMat[which(PiMat == 0)] <- 1e-6; logPi <- log(PiMat)
        AlphaNew <- UpdateBlockParameters(ExpectNew, dat$A)
        NewObjective <- ObjectiveFunction(dat$X, StoreMean, StoreCov, AlphaNew, logPi, dat$A, ExpectNew)
        if (NewObjective > Objective) {
          Alpha <- AlphaNew; Expect <- ExpectNew
          diff_val <- abs(Objective - NewObjective); Objective <- NewObjective
        } else { proceed <- 0 }
      }
      CommAssn <- apply(Expect, 1, function(x) which(x == max(x))[1])
      list(Comm = CommAssn, timed_out = FALSE, iters = iter, K = k)
    }, error = function(e) {
      if (grepl("reached elapsed time|reached CPU time", e$message)) {
        CommAssn <- apply(Expect, 1, function(x) which(x == max(x))[1])
        return(list(Comm = CommAssn, timed_out = TRUE, iters = iter, K = k))
      } else {
        return(list(Comm = NULL, timed_out = FALSE, iters = iter, K = k, error_msg = e$message))
      }
    })

    elapsed <- system.time({ force(fit) })[["elapsed"]]  # already timed above, just capture
    ari <- safe_ari(dat$z, fit$Comm)
    data.frame(scenario = task$scenario, sigma = task$sigma, r = task$ratio,
               rep = task$rep, N = task$N, K_true = task$K, K_est = K_est,
               K_final = fit$K, iters = fit$iters, final_ARI = ari,
               elapsed_sec = elapsed, timed_out = fit$timed_out,
               stringsAsFactors = FALSE)
  }

  tasks <- build_tasks()
  cat(sprintf("Pure-aSBM mode=%s: %d tasks, %d cores, %dh timeout\n",
              mode, length(tasks), n_cores, timeout_sec %/% 3600))

  cl <- makeCluster(min(n_cores, length(tasks)), type = "PSOCK")
  clusterEvalQ(cl, { library(igraph); source("FitAttribute.R") })
  clusterExport(cl, c("simulate_asbm_data", "safe_ari", "worker_fn"))

  res_list <- parLapply(cl, tasks, worker_fn)
  stopCluster(cl)

  res <- do.call(rbind, res_list)
  prefix <- if(mode == "A") paste0("pure_asbm_A_sigma", sigma)
  else paste0("pure_asbm_B_r", r_fixed)
  write.csv(res, paste0(prefix, "_results.csv"), row.names = FALSE)

  smry <- do.call(rbind, lapply(split(res, ~ scenario + r + sigma), function(sub) {
    data.frame(scenario = sub$scenario[1], sigma = sub$sigma[1], r = sub$r[1],
               reps = nrow(sub), ARI_mean = mean(sub$final_ARI, na.rm = TRUE),
               ARI_sd = sd(sub$final_ARI, na.rm = TRUE),
               elapsed_mean_sec = mean(sub$elapsed_sec),
               n_timeout = sum(sub$timed_out), stringsAsFactors = FALSE)
  }))
  write.csv(smry, paste0(prefix, "_point_summary.csv"), row.names = FALSE)

  cat(sprintf("\n=== Pure-aSBM %s Summary ===\n", mode))
  print(smry, row.names = FALSE)
  cat(sprintf("Saved: %s_results.csv, %s_point_summary.csv\n", prefix, prefix))
}

# ---- Main ----
args <- commandArgs(trailingOnly = TRUE)
mode <- if(length(args) >= 1) args[1] else "baselines"
submode <- if(length(args) >= 2) toupper(args[2]) else "A"
param <- if(length(args) >= 3) as.numeric(args[3]) else NULL

if (mode == "all") {
  cat("=== Running all experiments ===\n")
  run_baselines("A", sigma = if(is.null(param)) 3.0 else param)
  run_baselines("B")
  run_pure_asbm("A", sigma = if(is.null(param)) 3.0 else param)
  run_pure_asbm("B", r_fixed = if(is.null(param)) 4.0 else param)
  cat("\nAll experiments complete.\n")
} else if (mode == "baselines") {
  sigma_val <- if(!is.null(param)) param else 3.0
  run_baselines(submode, sigma = sigma_val)
} else if (mode == "pure") {
  if (submode == "A") {
    run_pure_asbm("A", sigma = if(!is.null(param)) param else 3.0)
  } else {
    run_pure_asbm("B", r_fixed = if(!is.null(param)) param else 4.0)
  }
} else {
  cat("Usage: Rscript run_simulations.R [baselines|pure|all] [A|B] [sigma|r]\n")
}
