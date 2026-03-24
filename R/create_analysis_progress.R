#' Create a Stage-Based Analysis Progress Reporter
#'
#' @param enabled Logical scalar indicating whether progress messages should be
#'   emitted.
#' @param total_steps Integer number of major workflow steps.
#'
#' @return A list containing a `tick()` function that advances the stage
#'   counter.
#' @keywords internal
#' @noRd
create_analysis_progress <- function(enabled, total_steps) {
  current_step <- 0L

  list(
    tick = function(label) {
      current_step <<- current_step + 1L

      if (isTRUE(enabled)) {
        message(sprintf("[%d/%d] %s", current_step, total_steps, label))
      }

      invisible(current_step)
    }
  )
}
