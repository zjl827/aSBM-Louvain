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
  Denom <- outer(N_k, N_k)
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
  LL_A_Val <- sum(q_curr * GraphLL_Matrix)

  log_q <- log(q_curr + 1e-20)
  log_pi <- log(PriorPi + 1e-20)
  KL_Val <- sum(q_curr * (log_q - log_pi))

  Total_ELBO <- LL_A_Val + LL_X_Val - Tau * KL_Val
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

# ------------------------------------------------------------------------------
# Main Algorithm: aSBM-Louvain
# ------------------------------------------------------------------------------

aSBM_Louvain <- function(AdjMat, AttributeMat,
                         T_outer = 20,
                         rho_tau = 0.8,
                         tau0 = 1.0,
                         resolution = 1.0,
                         epsilon_Q = 1e-4,
                         epsilon_L = 1e-3,
                         init_method = "kmeans"){

  N <- nrow(AttributeMat)

  # ---- Step 1: Initialization ----
  cat(sprintf("--- Initialization (%s) ---\n", init_method))

  if(is.null(colnames(AdjMat))) colnames(AdjMat) <- 1:N
  graph_obj <- graph_from_adjacency_matrix(AdjMat, mode = "undirected")

  if(init_method == "kmeans"){
    temp_louvain <- cluster_louvain(graph_obj, resolution = resolution)
    K_est <- length(unique(membership(temp_louvain)))
    if(K_est < 2) K_est <- 2
    km_res <- kmeans(AttributeMat, centers = K_est, nstart = 5)
    z_curr <- km_res$cluster
  } else {
    louvain_init <- cluster_louvain(graph_obj, resolution = resolution)
    z_curr <- membership(louvain_init)
  }

  z_curr <- as.numeric(factor(z_curr))
  K <- max(z_curr)

  StoreMean <- list()
  StoreCov <- list()
  for(k in 1:K){
    ids <- which(z_curr == k)
    if(length(ids) > 0){
      StoreMean[[k]] <- colMeans(AttributeMat[ids, , drop = FALSE])
      StoreCov[[k]] <- EstCovMat(AttributeMat[ids, , drop = FALSE])
    } else {
      StoreMean[[k]] <- rep(0, ncol(AttributeMat))
      StoreCov[[k]] <- diag(1, ncol(AttributeMat))
    }
  }

  q_curr <- matrix(0, nrow = N, ncol = K)
  for(i in 1:N) q_curr[i, z_curr[i]] <- 1

  p_curr <- UpdateBlockParameters(q_curr, AdjMat)
  pi_curr <- matrix(1 / K, nrow = N, ncol = K)

  tau <- tau0
  modularity_history <- c(modularity(
    graph_obj, membership = factor(z_curr),
    resolution = resolution, directed = FALSE
  ))
  elbo_history <- c(-Inf)

  # ---- Step 2: Main Loop ----
  for(t in 1:T_outer){

    # Phase A: Variational Update
    if(tau > 0){
      if(K == 1){
        pi_tilde <- matrix(1, nrow = N, ncol = 1)
      } else {
        pi_tilde <- matrix(0.01 / (K - 1), nrow = N, ncol = K)
        valid_z <- z_curr
        valid_z[valid_z > K] <- sample(1:K, sum(valid_z > K), replace = TRUE)
        for(i in 1:N) pi_tilde[i, valid_z[i]] <- 0.99
        pi_tilde <- pi_tilde / rowSums(pi_tilde)
      }
      pi_curr <- pi_tilde
    }

    q_curr <- VariationalEStep_Strict(AttributeMat, StoreMean, StoreCov,
                                       p_curr, pi_curr, tau, AdjMat, q_curr)
    q_curr <- q_curr / rowSums(q_curr)
    p_curr <- UpdateBlockParameters(q_curr, AdjMat)

    for(k in 1:K){
      StoreMean[[k]] <- MuUpdate(q_curr, AttributeMat, k)
      StoreCov[[k]] <- CovMatUpdate(StoreMean[[k]], q_curr, AttributeMat, k)
    }

    # Phase B: Louvain Refinement
    edges <- as_edgelist(graph_obj, names = FALSE)
    weights <- numeric(nrow(edges))
    QP <- q_curr %*% p_curr
    for(e_idx in 1:nrow(edges)){
      u <- edges[e_idx, 1]; v <- edges[e_idx, 2]
      weights[e_idx] <- sum(QP[u, ] * q_curr[v, ])
    }

    W_graph <- graph_from_adjacency_matrix(AdjMat, mode = "undirected", weighted = TRUE)
    E(W_graph)$weight <- weights

    louvain_res <- cluster_louvain(W_graph, resolution = resolution)
    z_new <- membership(louvain_res)
    mod_new <- modularity(
      W_graph, membership = z_new, weights = E(W_graph)$weight,
      resolution = resolution, directed = FALSE
    )

    # Phase C: Convergence Check
    elbo_curr <- Calculate_ELBO(AttributeMat, StoreMean, StoreCov,
                                 p_curr, pi_curr, q_curr, AdjMat, tau)
    elbo_prev <- tail(elbo_history, 1)
    elbo_change <- abs(elbo_curr - elbo_prev) / (abs(elbo_prev) + 1e-10)
    mod_diff <- abs(mod_new - tail(modularity_history, 1))

    if(t > 1 && mod_diff < epsilon_Q && elbo_change < epsilon_L){
      z_curr <- z_new; break
    }

    modularity_history <- c(modularity_history, mod_new)
    elbo_history <- c(elbo_history, elbo_curr)
    z_curr <- z_new

    # Handle K change after Louvain
    K_new <- length(unique(z_curr))
    if(K_new != K){
      K <- K_new; z_curr <- as.numeric(factor(z_curr))
      q_curr <- matrix(0, nrow = N, ncol = K)
      for(i in 1:N) q_curr[i, z_curr[i]] <- 1
      pi_curr <- matrix(1 / K, nrow = N, ncol = K)
      p_curr <- UpdateBlockParameters(q_curr, AdjMat)
      StoreMean <- list(); StoreCov <- list()
      for(k in 1:K){
        ids <- which(z_curr == k)
        if(length(ids) > 0){
          StoreMean[[k]] <- colMeans(AttributeMat[ids, , drop = FALSE])
          StoreCov[[k]] <- EstCovMat(AttributeMat[ids, , drop = FALSE])
        } else {
          StoreMean[[k]] <- rep(0, ncol(AttributeMat))
          StoreCov[[k]] <- diag(1, ncol(AttributeMat))
        }
      }
    }
    tau <- tau * rho_tau
  }

  return(list(
    partition = as.vector(z_curr),
    soft_labels = q_curr
  ))
}
