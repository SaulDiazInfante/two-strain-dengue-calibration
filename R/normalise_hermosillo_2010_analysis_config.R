#' Normalize a Hermosillo 2010 Analysis Config
#'
#' @param config A nested configuration list.
#'
#' @return A normalized configuration list ready for analysis execution.
#' @keywords internal
#' @noRd
normalise_hermosillo_2010_analysis_config <- function(config) {
  if (!is.list(config)) {
    stop("Analysis config must be a list.", call. = FALSE)
  }

  output_write <- config$outputs$write_outputs

  if (is.null(output_write)) {
    output_write <- config$outputs$write
  }

  optimisation_start <- as.numeric(config$optimisation$start)
  if (length(optimisation_start) != 2L) {
    stop("The optimisation start vector must have length 2.", call. = FALSE)
  }

  list(
    data = list(path = config$data$path),
    observed_rows = as.integer(resolve_sequence_spec(config$observed_rows)),
    grids = list(
      time = list(
        full = resolve_sequence_spec(config$grids$time$full),
        calibration = resolve_sequence_spec(config$grids$time$calibration),
        fit = resolve_sequence_spec(config$grids$time$fit, allow_null = TRUE)
      ),
      parameters = list(
        beta_h = resolve_sequence_spec(config$grids$parameters$beta_h),
        beta_m = resolve_sequence_spec(config$grids$parameters$beta_m),
        r01 = resolve_sequence_spec(config$grids$parameters$r01),
        r02 = resolve_sequence_spec(config$grids$parameters$r02),
        r0 = resolve_sequence_spec(config$grids$parameters$r0)
      )
    ),
    outputs = list(
      dir = if (is.null(config$outputs$dir)) "outputs" else config$outputs$dir,
      timestamp = config$outputs$timestamp,
      write_outputs = resolve_config_flag(
        if (is.null(output_write)) TRUE else output_write
      ),
      show_progress = resolve_config_flag(
        if (is.null(config$outputs$show_progress)) {
          "interactive"
        } else {
          config$outputs$show_progress
        }
      )
    ),
    optimisation = list(start = optimisation_start)
  )
}
