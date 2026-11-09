# Isolated, unchanged baseline environment; beta changes only the two likelihood expressions.
make_model <- function(code_dir, beta = NULL) {
  e <- new.env(parent = .GlobalEnv)
  sys.source(file.path(code_dir, "asbmlou.R"), envir = e)
  sys.source(file.path(code_dir, "fit_asbm_louvain.R"), envir = e)
  if (!is.null(beta)) {
    stopifnot(length(beta) == 1L, is.finite(beta), beta >= 0)
    e$attribute_beta <- beta
    e$Calculate_ELBO <- function(AttributeMat, StoreMean, StoreCov, Alpha, PriorPi, q_curr, Network, Tau) {
      LL_X_Val <- 0
      for (k in seq_len(ncol(q_curr))) {
        ll_vec <- GaussianLikelihood(AttributeMat, StoreMean[[k]], StoreCov[[k]])
        LL_X_Val <- LL_X_Val + sum(q_curr[, k] * ll_vec)
      }
      GraphLL_Matrix <- StructuralMessageFullBernoulli(Network, q_curr, Alpha)
      LL_A_Val <- 0.5 * sum(q_curr * GraphLL_Matrix)
      KL_Val <- sum(q_curr * (log(q_curr + 1e-20) - log(PriorPi + 1e-20)))
      LL_A_Val + attribute_beta * LL_X_Val - max(Tau, 1e-5) * KL_Val
    }
    e$VariationalEStep_Strict <- function(AttributeMat, StoreMean, StoreCov, Alpha, PriorPi, Tau, Network, LastQ) {
      AttributeLL <- matrix(0, nrow(AttributeMat), ncol(LastQ))
      for (k in seq_len(ncol(LastQ))) {
        AttributeLL[, k] <- GaussianLikelihood(AttributeMat, StoreMean[[k]], StoreCov[[k]])
      }
      GraphLL <- StructuralMessageFullBernoulli(Network, LastQ, Alpha)
      LogUnnormalized <- log(PriorPi + 1e-20) +
        (GraphLL + attribute_beta * AttributeLL) / max(Tau, 1e-5)
      t(apply(LogUnnormalized, 1, function(x) {
        ex <- exp(x - max(x)); ex / sum(ex)
      }))
    }
    environment(e$Calculate_ELBO) <- e
    environment(e$VariationalEStep_Strict) <- e
  }
  e
}
