#' Plot a Model Fit Against Observed Cases
#'
#' @description
#' Overlays simulated reported cases and observed Hermosillo case counts over a
#' selected time window.
#'
#' @param model_output Matrix or data frame returned by [run_model()].
#' @param data A data frame with columns `week`, `classical`, and
#'   `hemorrhagic`.
#' @param rows Integer row indices to display.
#' @param title Plot title.
#' @param week_origin Week number corresponding to simulation time zero.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#'
#' @return Invisibly returns the fitted case matrix used in the overlay.
#' @family hermosillo_analysis
#' @export
plot_model_fit <- function(model_output,
                           data,
                           rows = seq_len(nrow(data)),
                           title = "Model Fit",
                           week_origin = min(data$week[rows]),
                           xlim = NULL,
                           ylim = NULL) {
  fitted_cases <- extract_reported_cases(model_output)
  selected <- data[rows, , drop = FALSE]

  if (is.null(ylim)) {
    ylim <- c(0, max(fitted_cases, selected$classical, selected$hemorrhagic))
  }

  graphics::matplot(
    model_output[, "time"] + week_origin,
    fitted_cases,
    type = "l",
    lty = 1,
    lwd = 2,
    col = c(1, 2),
    cex.lab = 1.25,
    cex.axis = 1.25,
    main = title,
    xlab = "Week",
    ylab = "Infecteds",
    cex.main = 1.30,
    xlim = xlim,
    ylim = ylim
  )
  graphics::legend(
    "topleft",
    c("Classical", "Hemorrhagic"),
    lty = c(1, 1),
    col = c(1, 2),
    bty = "n",
    y.intersp = 0.35,
    inset = -0.025
  )
  graphics::points(selected$week, selected$classical, pch = 19, col = 1, cex = 0.8)
  graphics::points(selected$week, selected$hemorrhagic, pch = 19, col = 2, cex = 0.8)

  invisible(fitted_cases)
}
