#' Create a Text Progress Bar Descriptor
#'
#' @param enabled Logical scalar indicating whether a text progress bar should
#'   be created.
#' @param total Integer number of progress steps.
#' @param label Optional character label printed before the bar.
#'
#' @return A list describing the progress-bar state.
#' @keywords internal
#' @noRd
create_text_progress_bar <- function(enabled, total, label = NULL) {
  if (!isTRUE(enabled)) {
    return(list(enabled = FALSE, total = total, bar = NULL))
  }

  if (!is.null(label)) {
    message(label)
  }

  list(
    enabled = TRUE,
    total = total,
    bar = utils::txtProgressBar(min = 0, max = total, style = 3, file = stderr())
  )
}
