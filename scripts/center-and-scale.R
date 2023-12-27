library(purrr)
library(tibble)

center_and_scale_covs <- function(x_design_matrix, weights, scaleIndicators) {
  maybe_center_and_scale <- function(x, scaleIndicator) {
    `if`(
      scaleIndicator,
      center_and_scale(x),
      list(new_x = x, scale = 1)
    )
  }
  center_and_scale <- function(x) {
    mean_x <- sum(weights * x) / sum_weights
    centered_x <- x - mean_x
    weighted_abs_sum <- sum(abs(centered_x * weights))
    stopifnot(weighted_abs_sum > 0)
    scale_val <- sum_weights / weighted_abs_sum
    list(new_x = centered_x * scale_val, scale = scale_val)
  }
  stopifnot(
    is.data.frame(x_design_matrix),
    nrow(x_design_matrix) >= 1L,
    map_lgl(x_design_matrix, is.double),
    is.double(weights),
    length(weights) == nrow(x_design_matrix),
    map_lgl(scaleIndicators, is.logical),
    length(scaleIndicators) == length(x_design_matrix)
  )
  sum_weights <- sum(weights)
  stopifnot(sum_weights != 0)
  results <- map2(x_design_matrix, scaleIndicators, maybe_center_and_scale)
  list(
    x_design_matrix = as_tibble(map(results, "new_x")),
    scale_vals = map_dbl(results, "scale")
  )
}

covs <- tibble(
  x1 = c(1, 2, 3, 4, 5, 6),
  x2 = c(0, 1, 1, 0, 0, 2)
)
const_weights <- rep(1, nrow(covs))
nonconst_weights <- seq(1, 6, 1)


v1 <- center_and_scale_covs(covs, const_weights, c(FALSE, FALSE))

out <- center_and_scale_covs(covs, const_weights, c(TRUE, FALSE))
out_x1 <- out$x_design_matrix$x1
out_x2 <- out$x_design_matrix$x2
out_scale1 <- out$scale_vals["x1"]
sum_weights <- sum(const_weights)
centered <- covs$x1 - mean(covs$x1)
scale_val <- set_names(sum_weights / sum(abs(centered)), "x1")
expected <- centered * scale_val
all.equal(out_x1, expected)
all.equal(out_x2, covs$x2)

out <- center_and_scale_covs(covs, nonconst_weights, c(TRUE, FALSE))
out_x1 <- out$x_design_matrix$x1
out_scale1 <- out$scale_vals["x1"]
sum_weights <- sum(nonconst_weights)
weighted_mean <- sum(nonconst_weights * covs$x1) / sum_weights
centered <- covs$x1 - weighted_mean
scale_val <- set_names(sum_weights / sum(abs(nonconst_weights * centered)), "x1")
expected <- centered * scale_val
all.equal(out_x1, expected)
all.equal(out_scale1, scale_val)
