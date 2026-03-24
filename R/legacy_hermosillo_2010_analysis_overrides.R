#' Convert Legacy Hermosillo Runner Arguments into Config Overrides
#'
#' @param arguments Named list captured from `...`.
#'
#' @return A nested config override list.
#' @keywords internal
#' @noRd
legacy_hermosillo_2010_analysis_overrides <- function(arguments) {
  known_names <- c(
    "observed_rows",
    "full_time_grid",
    "calibration_grid",
    "fit_time_grid",
    "beta_h_values",
    "beta_m_values",
    "r01_values",
    "r02_values",
    "r0_values",
    "output_dir",
    "timestamp",
    "write_outputs",
    "show_progress",
    "optim_start"
  )

  if (is.null(names(arguments)) || any(names(arguments) == "")) {
    stop("Legacy overrides must be named arguments.", call. = FALSE)
  }

  unknown_names <- setdiff(names(arguments), known_names)
  if (length(unknown_names) > 0L) {
    stop(
      "Unsupported legacy arguments: ",
      paste(unknown_names, collapse = ", "),
      call. = FALSE
    )
  }

  drop_null_entries(list(
    observed_rows = arguments$observed_rows,
    grids = list(
      time = list(
        full = arguments$full_time_grid,
        calibration = arguments$calibration_grid,
        fit = arguments$fit_time_grid
      ),
      parameters = list(
        beta_h = arguments$beta_h_values,
        beta_m = arguments$beta_m_values,
        r01 = arguments$r01_values,
        r02 = arguments$r02_values,
        r0 = arguments$r0_values
      )
    ),
    outputs = list(
      dir = arguments$output_dir,
      timestamp = arguments$timestamp,
      write_outputs = arguments$write_outputs,
      show_progress = arguments$show_progress
    ),
    optimisation = list(start = arguments$optim_start)
  ))
}
