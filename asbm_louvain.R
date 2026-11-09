# ==============================================================================
# aSBM-Louvain: Attributed Network Clustering with Louvain-Based Block Structure Learning
# ==============================================================================
library('igraph')
library('mvtnorm')

# ------------------------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------------------------

GaussianLikelihood <- function(DM, mean, cov){
  mean <- as.vector(mean)
  tryCatch({
    dmvnorm(DM, mean = mean, sigma = cov, log = TRUE)
  }, error = function(e) {
    cov_fix <- cov + diag(1e-6, ncol(cov))
    dmvnorm(DM, mean = mean, sigma = cov_fix, log = TRUE)
  })
}

UpdateBlockParameters <- function(ExpMat, Network){
  Num <- t(ExpMat) %*% Network %*% ExpMat
  N_k <- colSums(ExpMat)
  Denom <- outer(N_k, N_k) - crossprod(ExpMat)
  NewBlockM <- Num / (Denom + 1e-10)
  NewBlockM <- (NewBlockM + t(NewBlockM)) / 2
  NewBlockM[NewBlockM > 0.999] <- 0.999
  NewBlockM[NewBlockM < 0.001] <- 0.001
  return(NewBlockM)
}

MuUpdate <- function(ExpMat, DataMat, commind){
  denom <- sum(ExpMat[, commind]) + 1e-10
  num <- colSums(DataMat * ExpMat[, commind])
  return(num / denom)
}

CovMatUpdate <- function(Mu, ExpMat, DataMat, commind){
  Mu <- as.vector(Mu)
  Diff <- t(t(DataMat) - Mu)
  w <- ExpMat[, commind]
  denom <- sum(w) + 1e-10
  CovMat <- (t(Diff) %*% (Diff * w)) / denom
  CovMat <- CovMat + diag(1e-6, ncol(DataMat))
  return(CovMat)
}

EstCovMat <- function(MVSamples){
  if(nrow(MVSamples) < 2) return(diag(1, ncol(MVSamples)))
  cov(MVSamples) + diag(1e-6, ncol(MVSamples))
}

# ------------------------------------------------------------------------------
# Full-Bernoulli Structural Message & ELBO
# ------------------------------------------------------------------------------

StructuralMessageFullBernoulli <- function(Network, LastQ, Alpha){
  LogP <- log(Alpha + 1e-10)
  EdgeMsg <- as.matrix(Network %*% (LastQ %*% LogP))

  LogOneMinusP <- log(1 - Alpha + 1e-10)
  NonedgeComponent <- LastQ %*% LogOneMinusP
  EdgeNonedgeComponent <- as.matrix(Network %*% NonedgeComponent)
  AllNodesComponent <- matrix(
    rep(colSums(NonedgeComponent), each = nrow(Network)),
    nrow = nrow(Network)
  )

  EdgeMsg + AllNodesComponent - NonedgeComponent - EdgeNonedgeComponent
}

Calculate_ELBO <- function(AttributeMat, StoreMean, StoreCov, Alpha, PriorPi, q_curr, Network, Tau){
  K <- ncol(q_curr)

  LL_X_Val <- 0
  for(k in 1:K){
    ll_vec <- GaussianLikelihood(AttributeMat, StoreMean[[k]], StoreCov[[k]])
    LL_X_Val <- LL_X_Val + sum(q_curr[, k] * ll_vec)
  }

  GraphLL_Matrix <- StructuralMessageFullBernoulli(Network, q_curr, Alpha)
  LL_A_Val <- 0.5 * sum(q_curr * GraphLL_Matrix)

  log_q <- log(q_curr + 1e-20)
  log_pi <- log(PriorPi + 1e-20)
  KL_Val <- sum(q_curr * (log_q - log_pi))

  Total_ELBO <- LL_A_Val + LL_X_Val - max(Tau, 1e-5) * KL_Val
  return(Total_ELBO)
}

# ------------------------------------------------------------------------------
# Variational E-Step
# ------------------------------------------------------------------------------

VariationalEStep_Strict <- function(AttributeMat, StoreMean, StoreCov, Alpha, PriorPi, Tau, Network, LastQ){
  N <- nrow(AttributeMat)
  K <- ncol(LastQ)

  AttributeLL <- matrix(0, nrow = N, ncol = K)
  for(k in 1:K){
    AttributeLL[, k] <- GaussianLikelihood(AttributeMat, StoreMean[[k]], StoreCov[[k]])
  }

  GraphLL <- StructuralMessageFullBernoulli(Network, LastQ, Alpha)

  SafeTau <- max(Tau, 1e-5)
  LogUnnormalized <- log(PriorPi + 1e-20) + (GraphLL + AttributeLL) / SafeTau

  NextQ <- t(apply(LogUnnormalized, 1, function(x) {
    mx <- max(x)
    ex <- exp(x - mx)
    return(ex / sum(ex))
  }))

  return(NextQ)
}

align_partition_max_overlap <- function(reference, candidate) {
  reference <- as.numeric(factor(reference))
  candidate <- as.numeric(factor(candidate))
  k_ref <- max(reference)
  k_cand <- max(candidate)

  if (k_ref != k_cand || k_ref <= 1L) {
    return(candidate)
  }

  overlap <- table(
    factor(reference, levels = seq_len(k_ref)),
    factor(candidate, levels = seq_len(k_cand))
  )

  edge_frame <- expand.grid(
    old = paste0("old_", seq_len(k_ref)),
    new = paste0("new_", seq_len(k_cand)),
    stringsAsFactors = FALSE
  )
  edge_frame$weight <- as.numeric(overlap[cbind(
    as.integer(sub("old_", "", edge_frame$old)),
    as.integer(sub("new_", "", edge_frame$new))
  )])

  alignment_graph <- igraph::graph_from_data_frame(edge_frame, directed = FALSE)
  vertex_types <- startsWith(igraph::V(alignment_graph)$name, "new_")
  matching <- igraph::max_bipartite_match(
    alignment_graph,
    types = vertex_types,
    weights = igraph::E(alignment_graph)$weight
  )$matching

  candidate_to_reference <- integer(k_cand)
  for (candidate_label in seq_len(k_cand)) {
    vertex_name <- paste0("new_", candidate_label)
    matched_vertex <- matching[vertex_name]
    if (is.na(matched_vertex)) {
      candidate_to_reference[candidate_label] <- candidate_label
    } else {
      candidate_to_reference[candidate_label] <- as.integer(
        sub("old_", "", unname(matched_vertex))
      )
    }
  }

  unname(candidate_to_reference[candidate])
}

aSBM_Louvain_revised <- function(
    AdjMat,
    AttributeMat,
    T_outer = 1000L,
    rho_tau = 0.8,
    tau0 = 1.0,
    resolution = 1.0,
    epsilon_Q = 1e-4,
    epsilon_L = 1e-3,
    k_stable_required = 3L,
    init_method = "kmeans",
    initial_partition = NULL,
    prior_feedback = TRUE,
    verbose = FALSE) {

  if (!init_method %in% c("kmeans", "louvain")) {
    stop("init_method must be 'kmeans' or 'louvain'.")
  }
  if (T_outer < 1L || k_stable_required < 1L) {
    stop("T_outer and k_stable_required must be positive integers.")
  }

  N <- nrow(AttributeMat)
  if (is.null(rownames(AdjMat))) rownames(AdjMat) <- seq_len(N)
  if (is.null(colnames(AdjMat))) colnames(AdjMat) <- rownames(AdjMat)
  graph_obj <- igraph::graph_from_adjacency_matrix(
    AdjMat, mode = "undirected", diag = FALSE
  )

  if (verbose) message(sprintf("Initialization: %s", init_method))
  if (!is.null(initial_partition)) {
    if (length(initial_partition) != N || anyNA(initial_partition)) {
      stop("initial_partition must contain one non-missing label per node.")
    }
    z_curr <- as.numeric(factor(initial_partition))
  } else if (init_method == "kmeans") {
    temp_louvain <- igraph::cluster_louvain(graph_obj, resolution = resolution)
    K_est <- length(unique(igraph::membership(temp_louvain)))
    if (K_est < 2L) K_est <- 2L
    z_curr <- stats::kmeans(AttributeMat, centers = K_est, nstart = 5)$cluster
  } else {
    z_curr <- igraph::membership(
      igraph::cluster_louvain(graph_obj, resolution = resolution)
    )
  }

  z_curr <- as.numeric(factor(z_curr))
  initial_partition_used <- z_curr
  K <- max(z_curr)

  initialise_parameters <- function(labels) {
    labels <- as.numeric(factor(labels))
    k_now <- max(labels)
    means <- vector("list", k_now)
    covariances <- vector("list", k_now)
    for (community in seq_len(k_now)) {
      ids <- which(labels == community)
      if (length(ids) > 0L) {
        means[[community]] <- colMeans(AttributeMat[ids, , drop = FALSE])
        covariances[[community]] <- EstCovMat(
          AttributeMat[ids, , drop = FALSE]
        )
      } else {
        means[[community]] <- rep(0, ncol(AttributeMat))
        covariances[[community]] <- diag(1, ncol(AttributeMat))
      }
    }
    q <- matrix(0, nrow = N, ncol = k_now)
    q[cbind(seq_len(N), labels)] <- 1
    list(
      labels = labels,
      K = k_now,
      means = means,
      covariances = covariances,
      q = q,
      p = UpdateBlockParameters(q, AdjMat),
      pi = matrix(1 / k_now, nrow = N, ncol = k_now)
    )
  }

  state <- initialise_parameters(z_curr)
  K <- state$K
  StoreMean <- state$means
  StoreCov <- state$covariances
  q_curr <- state$q
  p_curr <- state$p
  pi_curr <- state$pi

  tau <- tau0
  modularity_history <- numeric(0)
  elbo_history <- numeric(0)
  elbo_change_history <- numeric(0)
  modularity_change_history <- numeric(0)
  K_history <- K
  k_stable_run <- 0L
  converged <- FALSE
  iterations_run <- 0L
  stop_reason <- "safety_cap"

  edges <- igraph::as_edgelist(graph_obj, names = FALSE)

  for (t in seq_len(T_outer)) {
    iterations_run <- t

    if (!prior_feedback) {
      pi_curr <- matrix(1 / K, nrow = N, ncol = K)
    } else if (K == 1L) {
      pi_curr <- matrix(1, nrow = N, ncol = 1L)
    } else {
      pi_curr <- matrix(0.01 / (K - 1), nrow = N, ncol = K)
      valid_z <- z_curr
      valid_z[valid_z > K] <- sample(
        seq_len(K), sum(valid_z > K), replace = TRUE
      )
      pi_curr[cbind(seq_len(N), valid_z)] <- 0.99
      pi_curr <- pi_curr / rowSums(pi_curr)
    }

    q_curr <- VariationalEStep_Strict(
      AttributeMat, StoreMean, StoreCov, p_curr, pi_curr, tau, AdjMat, q_curr
    )
    q_curr <- q_curr / rowSums(q_curr)
    p_curr <- UpdateBlockParameters(q_curr, AdjMat)

    for (community in seq_len(K)) {
      StoreMean[[community]] <- MuUpdate(
        q_curr, AttributeMat, community
      )
      StoreCov[[community]] <- CovMatUpdate(
        StoreMean[[community]], q_curr, AttributeMat, community
      )
    }

    QP <- q_curr %*% p_curr
    weights <- numeric(nrow(edges))
    for (edge_index in seq_len(nrow(edges))) {
      u <- edges[edge_index, 1]
      v <- edges[edge_index, 2]
      weights[edge_index] <- sum(QP[u, ] * q_curr[v, ])
    }

    weighted_graph <- igraph::graph_from_adjacency_matrix(
      AdjMat, mode = "undirected", weighted = TRUE, diag = FALSE
    )
    igraph::E(weighted_graph)$weight <- weights
    louvain_result <- igraph::cluster_louvain(
      weighted_graph,
      weights = igraph::E(weighted_graph)$weight,
      resolution = resolution
    )
    z_new <- as.numeric(factor(igraph::membership(louvain_result)))
    K_new <- length(unique(z_new))
    if (K_new == K) {
      z_new <- align_partition_max_overlap(z_curr, z_new)
    }

    modularity_new <- igraph::modularity(
      weighted_graph,
      membership = z_new,
      weights = igraph::E(weighted_graph)$weight,
      resolution = resolution,
      directed = FALSE
    )
    elbo_new <- Calculate_ELBO(
      AttributeMat, StoreMean, StoreCov, p_curr, pi_curr,
      q_curr, AdjMat, tau
    )

    if (length(elbo_history) == 0L || K_new != tail(K_history, 1L)) {
      elbo_change <- Inf
    } else {
      elbo_previous <- tail(elbo_history, 1L)
      elbo_change <- abs(elbo_new - elbo_previous) /
        max(abs(elbo_previous), 1e-10)
    }
    modularity_change <- if (length(modularity_history) == 0L) {
      Inf
    } else {
      abs(modularity_new - tail(modularity_history, 1L))
    }

    if (K_new == tail(K_history, 1L)) {
      k_stable_run <- k_stable_run + 1L
    } else {
      k_stable_run <- 0L
    }

    elbo_history <- c(elbo_history, elbo_new)
    modularity_history <- c(modularity_history, modularity_new)
    elbo_change_history <- c(elbo_change_history, elbo_change)
    modularity_change_history <- c(
      modularity_change_history, modularity_change
    )
    K_history <- c(K_history, K_new)
    z_curr <- z_new

    if (
      is.finite(elbo_change) &&
      elbo_change < epsilon_L &&
      modularity_change < epsilon_Q &&
      k_stable_run >= k_stable_required
    ) {
      converged <- TRUE
      stop_reason <- "joint_stability"
      break
    }

    if (K_new != K) {
      state <- initialise_parameters(z_curr)
      K <- state$K
      StoreMean <- state$means
      StoreCov <- state$covariances
      q_curr <- state$q
      p_curr <- state$p
      pi_curr <- state$pi
    }

    tau <- tau * rho_tau
  }

  list(
    partition = as.vector(z_curr),
    initial_partition = as.vector(initial_partition_used),
    initial_K = length(unique(initial_partition_used)),
    soft_labels = q_curr,
    iterations_run = iterations_run,
    converged = converged,
    stop_reason = stop_reason,
    final_K = length(unique(z_curr)),
    prior_feedback = prior_feedback,
    final_elbo_change = tail(elbo_change_history, 1L),
    final_modularity_change = tail(modularity_change_history, 1L),
    final_k_stable_run = k_stable_run,
    elbo_history = elbo_history,
    modularity_history = modularity_history,
    K_history = K_history
  )
}
