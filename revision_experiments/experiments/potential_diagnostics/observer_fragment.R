    # Read-only diagnostics: no model assignments, stopping decisions or RNG draws.
    diagnostic_row <- data.frame(
      iteration=t, K_before=K, K_after=K_new, eligible=FALSE,
      reason='first_record', I=NA_real_, II=NA_real_, III=NA_real_, IV=NA_real_,
      delta_potential=NA_real_, epsilon=NA_real_, identity_residual=NA_real_,
      relative_identity_residual=NA_real_)
    if (!is.null(diagnostic_previous)) {
      if (K_new != K) {
        diagnostic_row$reason <- 'current_K_change'
      } else if (diagnostic_previous$K_before != diagnostic_previous$K_after ||
                 diagnostic_previous$K_after != K) {
        diagnostic_row$reason <- 'previous_K_change_or_dimension_mismatch'
      } else {
        old_context_value <- Calculate_ELBO(
          AttributeMat, StoreMean, StoreCov, p_curr,
          diagnostic_previous$prior, q_curr, AdjMat, diagnostic_previous$tau)
        old_partition_new_graph <- igraph::modularity(
          weighted_graph, membership=diagnostic_previous$partition,
          weights=igraph::E(weighted_graph)$weight, resolution=resolution,
          directed=FALSE)
        diagnostic_row$eligible <- TRUE
        diagnostic_row$reason <- 'fixed_dimension_aligned'
        diagnostic_row$I <- old_context_value-diagnostic_previous$objective
        diagnostic_row$II <- old_partition_new_graph-diagnostic_previous$modularity
        diagnostic_row$III <- modularity_new-old_partition_new_graph
        diagnostic_row$IV <- elbo_new-old_context_value
        diagnostic_row$delta_potential <-
          (elbo_new+diagnostic_alpha*modularity_new)-diagnostic_previous$potential
        parts <- c(diagnostic_row$I,diagnostic_alpha*diagnostic_row$II,
                   diagnostic_alpha*diagnostic_row$III,diagnostic_row$IV)
        diagnostic_row$epsilon <- sum(pmax(-parts,0))
        diagnostic_row$identity_residual <- diagnostic_row$delta_potential-sum(parts)
        diagnostic_row$relative_identity_residual <-
          abs(diagnostic_row$identity_residual)/max(1,sum(abs(parts)))
      }
    }
    diagnostic_rows[[length(diagnostic_rows)+1L]] <- diagnostic_row
    diagnostic_previous <- list(
      prior=pi_curr,tau=tau,partition=z_new,K_before=K,K_after=K_new,
      objective=elbo_new,modularity=modularity_new,
      potential=elbo_new+diagnostic_alpha*modularity_new)
