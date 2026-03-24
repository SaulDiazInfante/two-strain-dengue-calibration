#' Plot a Relative Likelihood Surface
#'
#' @description
#' Draws an image plot with a contour overlay for a relative likelihood surface
#' and optionally highlights a reference point.
#'
#' @param x_values Numeric vector for the x-axis.
#' @param y_values Numeric vector for the y-axis.
#' @param surface Numeric matrix of likelihood values.
#' @param x_label X-axis label.
#' @param y_label Y-axis label.
#' @param point Optional numeric vector of length 2 to highlight on the plot.
#' @param label Optional legend label for `point`.
#' @param main Plot title.
#'
#' @return Invisibly returns `surface`.
#' @family hermosillo_analysis
#' @export
plot_likelihood_surface <- function(x_values,
                                    y_values,
                                    surface,
                                    x_label,
                                    y_label,
                                    point = NULL,
                                    label = NULL,
                                    main = "Likelihood Contours") {
  fields::image.plot(
    x_values,
    y_values,
    surface,
    xlab = x_label,
    ylab = y_label,
    cex.lab = 1.25,
    main = main
  )
  graphics::contour(
    x_values,
    y_values,
    surface,
    levels = c(0.05),
    add = TRUE,
    drawlabels = TRUE,
    labcex = 1
  )

  if (!is.null(point)) {
    graphics::points(point[1], point[2], pch = 19, col = 2, cex = 1.4)
  }

  if (!is.null(label)) {
    graphics::legend(
      "topright",
      legend = label,
      pch = 19,
      col = 2,
      cex = 1.2,
      bty = "n",
      inset = -0.04
    )
  }

  invisible(surface)
}
