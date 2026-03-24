#' Format Profile-Likelihood Interval Labels
#'
#' @param interval Named numeric vector with entries `lower` and `upper`.
#'
#' @return A character vector of formatted labels, or `NULL` when the interval
#'   contains missing values.
#' @keywords internal
#' @noRd
format_interval_labels <- function(interval) {
  if (anyNA(interval)) {
    return(NULL)
  }

  c(
    sprintf("LI=%.4f", interval[["lower"]]),
    sprintf("LS=%.4f", interval[["upper"]])
  )
}
