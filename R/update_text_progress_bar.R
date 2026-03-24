#' Update a Text Progress Bar
#'
#' @param progress A progress descriptor created by
#'   `create_text_progress_bar()`.
#' @param value Integer progress value to display.
#'
#' @return Invisibly returns `progress`.
#' @keywords internal
#' @noRd
update_text_progress_bar <- function(progress, value) {
  if (isTRUE(progress$enabled)) {
    utils::setTxtProgressBar(progress$bar, value)
  }

  invisible(progress)
}
