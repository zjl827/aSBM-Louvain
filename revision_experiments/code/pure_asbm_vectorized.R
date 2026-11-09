suppressPackageStartupMessages(library(igraph))

pure_asbm_source_path <- tryCatch(
  normalizePath(sys.frame(1L)$ofile, mustWork = TRUE),
  error = function(condition) NA_character_
)
pure_asbm_native_path <- if (is.na(pure_asbm_source_path)) {
  NA_character_
} else {
  file.path(dirname(pure_asbm_source_path), "pure_asbm_floor_correction.so")
}
if (
  !is.na(pure_asbm_native_path) &&
    file.exists(pure_asbm_native_path) &&
    !is.loaded("pure_asbm_floor_correction_c")
) {
  dyn.load(pure_asbm_native_path)
}

pure_asbm_matched_kmeans_initialization <- function(
    network, attributes, resolution = 1.0, nstart = 5L) {
  graph_object <- igraph::graph_from_adjacency_matrix(
    network, mode = "undirected", diag = FALSE
  )
  louvain_partition <- igraph::membership(igraph::cluster_louvain(
    graph_object, resolution = resolution
  ))
  estimated_groups <- length(unique(louvain_partition))
  if (estimated_groups < 2L) estimated_groups <- 2L
  if (estimated_groups > nrow(attributes)) {
    stop("Louvain estimated more groups than available observations.")
  }
  kmeans_partition <- stats::kmeans(
    attributes,
    centers = estimated_groups,
    nstart = nstart
  )$cluster
  list(
    partition = as.numeric(factor(kmeans_partition)),
    estimated_groups = estimated_groups,
    louvain_partition = as.numeric(factor(louvain_partition)),
    resolution = resolution,
    nstart = nstart
  )
}

pure_asbm_row_logsumexp <- function(values) {
  maxima <- apply(values, 1L, max)
  maxima + log(rowSums(exp(values - maxima)))
}

pure_asbm_gaussian_loglik <- function(data, mean_vector, covariance) {
  mvtnorm::dmvnorm(
    data,
    mean = as.numeric(mean_vector),
    sigma = covariance,
    log = TRUE
  )
}

pure_asbm_graph_loglik <- function(network, expectation, alpha) {
  neighbor_mass <- network %*% expectation
  total_mass_except_self <- matrix(
    colSums(expectation),
    nrow = nrow(expectation),
    ncol = ncol(expectation),
    byrow = TRUE
  ) - expectation
  nonneighbor_mass <- total_mass_except_self - neighbor_mass
  neighbor_mass %*% t(log(alpha)) +
    nonneighbor_mass %*% t(log1p(-alpha))
}

pure_asbm_floor_correction <- function(
    left_weights, right_weights, floor_value = 1e-8) {
  # Compute
  #   sum_{i < j} [max(floor_value, left_i * right_j)
  #                - left_i * right_j]
  # exactly up to floating-point summation order. Scanning i from right to
  # left, Fenwick trees store the count and sum of right_j values for j > i.
  if (is.loaded("pure_asbm_floor_correction_c")) {
    return(.Call(
      "pure_asbm_floor_correction_c",
      as.double(left_weights),
      as.double(right_weights),
      as.double(floor_value)
    ))
  }

  number_of_nodes <- length(left_weights)
  if (number_of_nodes < 2L) return(0)

  sorted_values <- sort(unique(right_weights))
  value_ranks <- match(right_weights, sorted_values)
  tree_counts <- integer(length(sorted_values))
  tree_sums <- numeric(length(sorted_values))
  correction <- 0

  for (node_index in number_of_nodes:1L) {
    left_value <- left_weights[node_index]
    threshold <- if (left_value > 0) {
      floor_value / left_value
    } else {
      Inf
    }
    low_rank <- findInterval(threshold, sorted_values)
    if (low_rank > 0L) {
      low_count <- 0L
      low_sum <- 0
      query_index <- low_rank
      while (query_index > 0L) {
        low_count <- low_count + tree_counts[query_index]
        low_sum <- low_sum + tree_sums[query_index]
        query_index <- query_index - bitwAnd(query_index, -query_index)
      }
      correction <- correction +
        floor_value * low_count - left_value * low_sum
    }

    update_index <- value_ranks[node_index]
    while (update_index <= length(sorted_values)) {
      tree_counts[update_index] <- tree_counts[update_index] + 1L
      tree_sums[update_index] <-
        tree_sums[update_index] + right_weights[node_index]
      update_index <- update_index + bitwAnd(update_index, -update_index)
    }
  }
  correction
}

pure_asbm_update_block_parameters <- function(expectation, network) {
  # The reference code sums over i < j for q <= l, applies the lower bound
  # max(1e-8, tau_iq * tau_jl) to every denominator term, and mirrors the
  # result. Matrix products compute the unfloored sums, and a Fenwick-tree
  # correction restores the per-pair lower bound without materializing an
  # N-by-N matrix for each community pair.
  number_of_groups <- ncol(expectation)
  upper_network <- Matrix::triu(Matrix::Matrix(network, sparse = TRUE), k = 1L)
  numerators <- crossprod(expectation, upper_network %*% expectation)
  suffix_mass <- apply(expectation, 2L, function(values) {
    rev(cumsum(rev(values))) - values
  })
  denominators <- crossprod(expectation, suffix_mass)
  output <- matrix(0, nrow = number_of_groups, ncol = number_of_groups)
  for (q in seq_len(number_of_groups)) {
    for (l in q:number_of_groups) {
      denominator <- denominators[q, l]
      if (
        min(expectation[, q]) * min(expectation[, l]) < 1e-8
      ) {
        denominator <- denominator + pure_asbm_floor_correction(
          expectation[, q], expectation[, l], floor_value = 1e-8
        )
      }
      value <- numerators[q, l] / denominator
      output[q, l] <- value
      output[l, q] <- value
    }
  }
  output
}

pure_asbm_update_gaussian <- function(data, expectation, group_index) {
  weights <- expectation[, group_index]
  denominator <- sum(weights)
  mean_vector <- colSums(data * weights) / denominator
  centered <- sweep(data, 2L, mean_vector, FUN = "-")
  covariance <- crossprod(centered, centered * weights) / denominator
  list(mean = mean_vector, covariance = covariance)
}

pure_asbm_objective <- function(
    data, means, covariances, alpha, log_pi, network, expectation) {
  number_of_groups <- ncol(expectation)
  attribute_loglik <- vapply(
    seq_len(number_of_groups),
    function(group_index) {
      pure_asbm_gaussian_loglik(
        data, means[[group_index]], covariances[[group_index]]
      )
    },
    numeric(nrow(data))
  )
  attribute_loglik[attribute_loglik == 0] <- 1e-7
  graph_loglik <- pure_asbm_graph_loglik(network, expectation, alpha)
  sum(expectation * (
    attribute_loglik + graph_loglik +
      matrix(log_pi, nrow(data), number_of_groups, byrow = TRUE)
  ))
}

pure_asbm_expectation <- function(
    data, means, covariances, alpha, log_pi, old_expectation, network) {
  number_of_groups <- ncol(old_expectation)
  attribute_loglik <- vapply(
    seq_len(number_of_groups),
    function(group_index) {
      pure_asbm_gaussian_loglik(
        data, means[[group_index]], covariances[[group_index]]
      )
    },
    numeric(nrow(data))
  )
  attribute_loglik[attribute_loglik == 0] <- 1e-8
  graph_loglik <- pure_asbm_graph_loglik(
    network, old_expectation, alpha
  ) + matrix(log_pi, nrow(data), number_of_groups, byrow = TRUE)
  unnormalized <- attribute_loglik + graph_loglik
  output <- exp(unnormalized - pure_asbm_row_logsumexp(unnormalized))
  output[output == 0] <- exp(-10)
  output
}

fit_pure_asbm_vectorized <- function(
    network, attributes, initial_partition = NULL,
    tolerance = 1, timeout_sec = 3600, max_iterations = 1000L) {
  if (is.null(initial_partition)) {
    graph_object <- igraph::graph_from_adjacency_matrix(
      network, mode = "undirected"
    )
    initial_partition <- igraph::membership(
      igraph::cluster_louvain(graph_object)
    )
  }
  initial_partition <- as.numeric(factor(initial_partition))
  number_of_nodes <- nrow(attributes)
  number_of_groups <- max(initial_partition)

  means <- vector("list", number_of_groups)
  covariances <- vector("list", number_of_groups)
  for (group_index in seq_len(number_of_groups)) {
    ids <- which(initial_partition == group_index)
    means[[group_index]] <- colMeans(attributes[ids, , drop = FALSE])
    centered <- sweep(
      attributes[ids, , drop = FALSE], 2L, means[[group_index]], FUN = "-"
    )
    covariances[[group_index]] <- crossprod(centered) / length(ids)
    if (sum(covariances[[group_index]]) == 0) {
      covariances[[group_index]] <- covariances[[group_index]] + 1e-6
    }
  }

  expectation <- matrix(0.01, nrow = number_of_nodes, ncol = number_of_groups)
  expectation[cbind(seq_len(number_of_nodes), initial_partition)] <- 0.99
  alpha <- pure_asbm_update_block_parameters(expectation, network)
  log_pi <- log(colSums(expectation) / number_of_nodes)
  objective <- pure_asbm_objective(
    attributes, means, covariances, alpha, log_pi, network, expectation
  )

  iteration <- 0L
  difference <- Inf
  stop_reason <- "objective_not_increased"
  start_time <- proc.time()[["elapsed"]]
  repeat {
    if (difference <= tolerance) {
      stop_reason <- "tolerance"
      break
    }
    if (iteration >= max_iterations) {
      stop_reason <- "iteration_cap"
      break
    }
    if (proc.time()[["elapsed"]] - start_time > timeout_sec) {
      stop_reason <- "timeout"
      break
    }
    iteration <- iteration + 1L

    for (group_index in seq_len(number_of_groups)) {
      update <- pure_asbm_update_gaussian(
        attributes, expectation, group_index
      )
      means[[group_index]] <- update$mean
      covariances[[group_index]] <- update$covariance
    }
    new_expectation <- pure_asbm_expectation(
      attributes, means, covariances, alpha, log_pi, expectation, network
    )
    new_pi <- colSums(new_expectation) / number_of_nodes
    new_pi[new_pi == 0] <- 1e-6
    log_pi <- log(new_pi)
    new_alpha <- pure_asbm_update_block_parameters(new_expectation, network)
    new_objective <- pure_asbm_objective(
      attributes, means, covariances, new_alpha, log_pi,
      network, new_expectation
    )

    if (!is.finite(new_objective) || new_objective <= objective) {
      stop_reason <- "objective_not_increased"
      break
    }
    alpha <- new_alpha
    expectation <- new_expectation
    difference <- abs(objective - new_objective)
    objective <- new_objective
  }

  partition <- max.col(expectation, ties.method = "first")
  list(
    partition = partition,
    expectation = expectation,
    alpha = alpha,
    means = means,
    covariances = covariances,
    iterations = iteration,
    final_objective = objective,
    final_difference = difference,
    stop_reason = stop_reason,
    initial_partition = initial_partition,
    final_K = length(unique(partition))
  )
}
