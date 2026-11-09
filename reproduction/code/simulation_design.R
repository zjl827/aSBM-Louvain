# Data generators used for the balanced and unequal-community simulations.

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

calibrate_probabilities <- function(block_sizes, ratio, target_degree) {
  number_of_nodes <- sum(block_sizes)
  within_ordered_pairs <- sum(block_sizes * (block_sizes - 1))
  cross_ordered_pairs <- number_of_nodes^2 - sum(block_sizes^2)
  denominator <-
    ratio * within_ordered_pairs + cross_ordered_pairs
  p_out <- number_of_nodes * target_degree / denominator
  p_in <- ratio * p_out
  if (p_in > 1 || p_out > 1) {
    stop("Calibrated SBM probability exceeds one.")
  }
  expected_degree_by_block <-
    (block_sizes - 1) * p_in + (number_of_nodes - block_sizes) * p_out
  list(
    p_in = p_in,
    p_out = p_out,
    expected_degree_overall = sum(block_sizes * expected_degree_by_block) /
      number_of_nodes,
    expected_degree_by_block = expected_degree_by_block
  )
}

simulate_unbalanced_data <- function(
    block_sizes, D, ratio, sigma, target_degree, seed) {
  set.seed(seed)
  number_of_nodes <- sum(block_sizes)
  number_of_blocks <- length(block_sizes)
  probabilities <- calibrate_probabilities(
    block_sizes, ratio, target_degree
  )
  preference <- matrix(
    probabilities$p_out,
    nrow = number_of_blocks,
    ncol = number_of_blocks
  )
  diag(preference) <- probabilities$p_in
  graph_object <- igraph::sample_sbm(
    number_of_nodes,
    pref.matrix = preference,
    block.sizes = block_sizes,
    directed = FALSE,
    loops = FALSE
  )
  adjacency <- as.matrix(
    igraph::as_adjacency_matrix(graph_object, sparse = FALSE)
  )
  truth <- rep(seq_len(number_of_blocks), times = block_sizes)
  centers <- matrix(
    stats::rnorm(number_of_blocks * D),
    nrow = number_of_blocks,
    ncol = D
  )
  attributes <- matrix(0, nrow = number_of_nodes, ncol = D)
  for (community in seq_len(number_of_blocks)) {
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
  list(
    A = adjacency,
    X = attributes,
    z = truth,
    graph = graph_object,
    probabilities = probabilities
  )
}
