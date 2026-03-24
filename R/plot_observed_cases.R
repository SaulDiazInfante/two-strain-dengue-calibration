#' Plot Observed Hermosillo Case Counts
#'
#' @description
#' Creates a scatter plot of classical and hemorrhagic dengue cases for a
#' selected subset of the Hermosillo study data.
#'
#' @param data A data frame with columns `week`, `classical`, and
#'   `hemorrhagic`.
#' @param rows Integer row indices to display.
#' @param title Plot title.
#' @param ylim Optional y-axis limits.
#'
#' @return Invisibly returns the plotted subset of `data`.
#' @family hermosillo_analysis
#' @export
plot_observed_cases <- function(data,
                                rows = seq_len(nrow(data)),
                                title = "Observed Cases",
                                ylim = NULL) {
  selected <- data[rows, , drop = FALSE]

  if (is.null(ylim)) {
    ylim <- c(0, max(selected$classical, selected$hemorrhagic))
  }

  graphics::matplot(
    selected$week,
    cbind(selected$classical, selected$hemorrhagic),
    type = "p",
    pch = c(19, 19),
    col = c(1, 2),
    cex.lab = 1.35,
    cex.axis = 1.35,
    main = title,
    xlab = "Week",
    ylab = "Infecteds",
    cex.main = 1.30,
    ylim = ylim
  )
  graphics::legend(
    "topleft",
    c("Classical", "Hemorrhagic"),
    cex = 1.2,
    pch = c(19, 19),
    col = c(1, 2),
    bty = "n",
    y.intersp = 0.35,
    inset = -0.025
  )

  invisible(selected)
}
