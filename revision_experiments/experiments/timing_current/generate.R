simulate_asbm_data_revised <- function(
    N, K, D, ratio, sigma, target_degree, seed) {
  set.seed(seed)
  block_sizes <- rep(N %/% K, K)
  if (N %% K > 0L) {
    block_sizes[seq_len(N %% K)] <-
      block_sizes[seq_len(N %% K)] + 1L
  }
  group_size <- N / K
  p_out <- target_degree /
    (((group_size - 1) * ratio) + (N - group_size))
  p_in <- ratio * p_out
  preference <- matrix(p_out, K, K)
  diag(preference) <- p_in
  graph_obj <- igraph::sample_sbm(
    N,
    pref.matrix = preference,
    block.sizes = block_sizes,
    directed = FALSE,
    loops = FALSE
  )
  adjacency <- as.matrix(
    igraph::as_adjacency_matrix(graph_obj, sparse = FALSE)
  )
  truth <- rep(seq_len(K), times = block_sizes)
  centers <- matrix(stats::rnorm(K * D), nrow = K, ncol = D)
  attributes <- matrix(0, nrow = N, ncol = D)
  for (community in seq_len(K)) {
    ids <- which(truth == community)
    attributes[ids, ] <- matrix(
      stats::rnorm(
        length(ids) * D,
        mean = rep(centers[community, ], each = length(ids)),
        sd = sigma
      ),
      nrow = length(ids),
      byrow = FALSE
    )
  }
  list(A = adjacency, X = attributes, z = truth, graph = graph_obj)
}
