#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>

static int compare_double(const void *left, const void *right) {
  const double a = *(const double *) left;
  const double b = *(const double *) right;
  return (a > b) - (a < b);
}

static R_xlen_t lower_bound_double(
    const double *values, R_xlen_t length, double target) {
  R_xlen_t low = 0;
  R_xlen_t high = length;
  while (low < high) {
    const R_xlen_t middle = low + (high - low) / 2;
    if (values[middle] < target) {
      low = middle + 1;
    } else {
      high = middle;
    }
  }
  return low;
}

static R_xlen_t upper_bound_double(
    const double *values, R_xlen_t length, double target) {
  R_xlen_t low = 0;
  R_xlen_t high = length;
  while (low < high) {
    const R_xlen_t middle = low + (high - low) / 2;
    if (values[middle] <= target) {
      low = middle + 1;
    } else {
      high = middle;
    }
  }
  return low;
}

SEXP pure_asbm_floor_correction_c(
    SEXP left_sexp, SEXP right_sexp, SEXP floor_sexp) {
  if (TYPEOF(left_sexp) != REALSXP || TYPEOF(right_sexp) != REALSXP) {
    error("left_weights and right_weights must be double vectors");
  }
  const R_xlen_t number_of_nodes = XLENGTH(left_sexp);
  if (XLENGTH(right_sexp) != number_of_nodes) {
    error("left_weights and right_weights must have equal length");
  }
  if (number_of_nodes < 2) return ScalarReal(0.0);

  const double floor_value = asReal(floor_sexp);
  const double *left = REAL(left_sexp);
  const double *right = REAL(right_sexp);

  double *sorted = (double *) R_alloc(number_of_nodes, sizeof(double));
  for (R_xlen_t i = 0; i < number_of_nodes; ++i) {
    if (!R_FINITE(left[i]) || !R_FINITE(right[i]) ||
        left[i] < 0.0 || right[i] < 0.0) {
      error("weights must be finite and nonnegative");
    }
    sorted[i] = right[i];
  }
  qsort(sorted, number_of_nodes, sizeof(double), compare_double);

  R_xlen_t unique_length = 1;
  for (R_xlen_t i = 1; i < number_of_nodes; ++i) {
    if (sorted[i] != sorted[unique_length - 1]) {
      sorted[unique_length++] = sorted[i];
    }
  }

  double *tree_counts = (double *) R_alloc(unique_length + 1, sizeof(double));
  double *tree_sums = (double *) R_alloc(unique_length + 1, sizeof(double));
  for (R_xlen_t i = 0; i <= unique_length; ++i) {
    tree_counts[i] = 0.0;
    tree_sums[i] = 0.0;
  }

  double correction = 0.0;
  for (R_xlen_t reverse_index = number_of_nodes; reverse_index > 0;
       --reverse_index) {
    const R_xlen_t node_index = reverse_index - 1;
    const double left_value = left[node_index];
    const double threshold = left_value > 0.0
      ? floor_value / left_value
      : R_PosInf;
    R_xlen_t query_index = upper_bound_double(
      sorted, unique_length, threshold
    );
    double low_count = 0.0;
    double low_sum = 0.0;
    while (query_index > 0) {
      low_count += tree_counts[query_index];
      low_sum += tree_sums[query_index];
      query_index -= query_index & -query_index;
    }
    correction += floor_value * low_count - left_value * low_sum;

    R_xlen_t update_index = lower_bound_double(
      sorted, unique_length, right[node_index]
    ) + 1;
    while (update_index <= unique_length) {
      tree_counts[update_index] += 1.0;
      tree_sums[update_index] += right[node_index];
      update_index += update_index & -update_index;
    }
  }

  return ScalarReal(correction);
}
