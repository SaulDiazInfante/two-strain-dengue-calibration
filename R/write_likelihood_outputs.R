#' Write Likelihood Surfaces to CSV Files
#'
#' @param loglikelihood Numeric vector or matrix of log-likelihood values.
#' @param log_file Destination path for the raw log-likelihood values.
#' @param relative_file Destination path for the relative likelihood values.
#'
#' @return A numeric vector or matrix of relative likelihood values.
#' @keywords internal
#' @noRd
write_likelihood_outputs <- function(loglikelihood, log_file, relative_file) {
  relative <- relative_likelihood(loglikelihood)
  utils::write.csv(loglikelihood, log_file)
  utils::write.csv(relative, relative_file)
  relative
}
