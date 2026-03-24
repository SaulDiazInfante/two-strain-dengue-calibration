#' Evaluate a Likelihood Grid with Optional Progress Updates
#'
#' @param x_values Numeric vector for the first parameter axis.
#' @param y_values Numeric vector for the second parameter axis.
#' @param evaluator Function taking a two-element numeric vector and returning a
#'   negative log-likelihood.
#' @param show_progress Logical scalar indicating whether a text progress bar
#'   should be shown.
#' @param progress_label Optional character label printed before the bar.
#'
#' @return A matrix of log-likelihood values.
#' @keywords internal
#' @noRd
evaluate_loglikelihood_grid_with_progress <- function(x_values,
                                                     y_values,
                                                     evaluator,
                                                     show_progress = FALSE,
                                                     progress_label = NULL) {
  if (length(x_values) == 0L || length(y_values) == 0L) {
    return(matrix(numeric(0), nrow = length(x_values), ncol = length(y_values)))
  }

  progress <- create_text_progress_bar(
    enabled = show_progress,
    total = length(x_values),
    label = progress_label
  )
  on.exit(close_text_progress_bar(progress), add = TRUE)

  row_results <- lapply(seq_along(x_values), function(i) {
    row_values <- vapply(
      y_values,
      function(y_value) -evaluator(c(x_values[i], y_value)),
      numeric(1)
    )

    update_text_progress_bar(progress, i)
    row_values
  })

  matrix(
    unlist(row_results, use.names = FALSE),
    nrow = length(x_values),
    ncol = length(y_values),
    byrow = TRUE
  )
}
