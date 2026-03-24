#' Save a Plot to a JPEG File
#'
#' @param path Destination file path.
#' @param plot_call A function with no arguments that draws the plot.
#'
#' @return Invisibly returns `path`.
#' @keywords internal
#' @noRd
save_plot <- function(path, plot_call) {
  grDevices::jpeg(path, width = 1800, height = 1200, res = 200)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_call()
  invisible(path)
}
