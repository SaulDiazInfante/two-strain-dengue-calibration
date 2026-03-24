#' Close a Text Progress Bar
#'
#' @param progress A progress descriptor created by
#'   `create_text_progress_bar()`.
#'
#' @return Invisibly returns `NULL`.
#' @keywords internal
#' @noRd
close_text_progress_bar <- function(progress) {
  if (isTRUE(progress$enabled)) {
    close(progress$bar)
  }

  invisible(NULL)
}
