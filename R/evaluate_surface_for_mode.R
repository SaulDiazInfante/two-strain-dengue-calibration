#' Evaluate a Reparameterised Likelihood Surface
#'
#' @param beta_values Numeric vector of `betaH` values.
#' @param reproduction_values Numeric vector of reproduction-number values.
#' @param mode Character scalar selecting `"R01"`, `"R02"`, or `"R0"`.
#' @param observation_times Observation times on the simulation grid.
#' @param observed_cases Two-column matrix of observed classical and
#'   hemorrhagic case counts.
#' @param initial_state Named numeric state vector used by [run_model()].
#' @param time_grid Numeric vector of solver times.
#' @param fixed_parameters Named numeric vector containing the fixed model
#'   parameters.
#' @param show_progress Logical scalar indicating whether a text progress bar
#'   should be shown.
#' @param progress_label Optional character label printed before the bar.
#' @param log_file Optional output path for the raw log-likelihood surface.
#' @param relative_file Optional output path for the relative likelihood
#'   surface.
#'
#' @return A list with the raw surface, relative surface, and written file
#'   paths.
#' @keywords internal
#' @noRd
evaluate_surface_for_mode <- function(beta_values,
                                      reproduction_values,
                                      mode,
                                      observation_times,
                                      observed_cases,
                                      initial_state,
                                      time_grid,
                                      fixed_parameters,
                                      show_progress = FALSE,
                                      progress_label = NULL,
                                      log_file = NULL,
                                      relative_file = NULL) {
  likelihood_fn <- switch(
    mode,
    R01 = negative_loglikelihood_r01,
    R02 = negative_loglikelihood_r02,
    R0 = negative_loglikelihood_r0,
    stop("Unsupported reproduction mode: ", mode, call. = FALSE)
  )

  loglikelihood <- evaluate_loglikelihood_grid_with_progress(
    beta_values,
    reproduction_values,
    function(par) {
      likelihood_fn(
        par = par,
        observation_times = observation_times,
        observed_cases = observed_cases,
        initial_state = initial_state,
        time_grid = time_grid,
        fixed_parameters = fixed_parameters
      )
    },
    show_progress = show_progress,
    progress_label = progress_label
  )

  relative <- if (is.null(log_file) || is.null(relative_file)) {
    relative_likelihood(loglikelihood)
  } else {
    write_likelihood_outputs(loglikelihood, log_file, relative_file)
  }

  list(
    loglikelihood = loglikelihood,
    relative = relative,
    files = list(loglikelihood = log_file, relative = relative_file)
  )
}
