# Reconstruct the first predetermined replication of each W01--W08
# Gaussian condition, with the same graph, means and residual streams.
gb_simulated_source <- function(source) {
  stopifnot(source$r %in% c(1.5, 4), source$N == 1000L, source$D == 50L)
  n <- 1000L; k <- 5L; d <- 50L; sizes <- rep(200L, k)
  truth <- rep(seq_len(k), times = sizes)
  ri <- match(source$r, c(1.5, 4))
  set.seed(92210000L + ri * 100000L + 1L)
  pout <- 20 / (199 * source$r + 800)
  block <- matrix(pout, k, k)
  diag(block) <- source$r * pout
  graph <- igraph::sample_sbm(n, block, sizes, directed = FALSE, loops = FALSE)
  A <- as.matrix(igraph::as_adjacency_matrix(graph, sparse = FALSE))
  set.seed(92211001L)
  means <- matrix(rnorm(k * d), k, d)
  set.seed(92212001L)
  Z <- matrix(rnorm(n * d), n, d)
  set.seed(92213001L)
  high <- sort(sample(seq_len(k), 2L))
  scales <- rep(1, k)
  scales[high] <- source$h
  sigmas <- source$sigma * scales / sqrt(sum(sizes / sum(sizes) * scales^2))
  X <- means[truth, , drop = FALSE] + Z * sigmas[truth]
  list(A = A, X = X, truth = truth)
}
