community_noise_sd <- function(h, sigma0, high_groups, block_sizes) {
  stopifnot(is.finite(h), h >= 1, sigma0 > 0,
            length(high_groups) == 2L, !anyDuplicated(high_groups),
            all(high_groups %in% seq_along(block_sizes)), all(block_sizes > 0))
  a <- rep(1, length(block_sizes))
  a[high_groups] <- h
  sigma0 * a / sqrt(sum(block_sizes / sum(block_sizes) * a^2))
}

partition_nmi <- function(truth, prediction) {
  tab <- table(truth, prediction)
  p <- tab / sum(tab)
  px <- rowSums(p); py <- colSums(p)
  expected <- outer(px, py)
  valid <- p > 0
  mi <- sum(p[valid] * log(p[valid] / expected[valid]))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx + hy == 0) return(1)
  max(0, min(1, 2 * mi / (hx + hy)))
}

community_scores <- function(truth, prediction, high_groups) {
  truth <- as.integer(factor(truth))
  prediction <- as.integer(factor(prediction))
  kt <- max(truth); kp <- max(prediction); km <- max(kt, kp)
  counts <- matrix(0, km, km)
  counts[seq_len(kt), seq_len(kp)] <- table(
    factor(truth, levels = seq_len(kt)),
    factor(prediction, levels = seq_len(kp))
  )
  edges <- expand.grid(old = paste0("t", seq_len(km)),
                       new = paste0("p", seq_len(km)),
                       stringsAsFactors = FALSE)
  edges$weight <- as.numeric(counts)
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  match <- igraph::max_bipartite_match(
    g, types = startsWith(igraph::V(g)$name, "p"),
    weights = igraph::E(g)$weight
  )$matching
  rows <- lapply(seq_len(kt), function(k) {
    matched <- unname(match[paste0("t", k)])
    j <- if (is.na(matched)) NA_integer_ else as.integer(sub("p", "", matched))
    true_size <- sum(truth == k)
    pred_size <- if (is.na(j) || j > kp) 0L else sum(prediction == j)
    tp <- if (is.na(j) || j > kp) 0L else counts[k, j]
    precision <- if (pred_size == 0) 0 else tp / pred_size
    recall <- tp / true_size
    f1 <- if (precision + recall == 0) 0 else 2 * precision * recall / (precision + recall)
    data.frame(community_id = k, precision = precision, recall = recall,
               F1 = f1, is_high_noise = k %in% high_groups,
               true_size = true_size, predicted_size = pred_size,
               matched_label = if (is.na(j) || j > kp) NA_integer_ else j)
  })
  do.call(rbind, rows)
}

hash_object <- function(x) {
  p <- tempfile("pilot_hash_")
  on.exit(unlink(p))
  saveRDS(x, p, compress = FALSE, version = 3)
  unname(tools::md5sum(p))
}
