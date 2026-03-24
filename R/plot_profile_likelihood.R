#' Plot a Relative Profile Likelihood
#'
#' @description
#' Plots a one-dimensional relative profile likelihood curve and optionally
#' annotates the corresponding confidence interval.
#'
#' @param values Numeric vector defining the profile grid.
#' @param relative_profile Numeric vector of relative profile likelihood values.
#' @param x_label X-axis label.
#' @param interval Named numeric vector with `lower` and `upper`.
#' @param interval_label Optional character vector for the interval legend.
#' @param main Plot title.
#'
#' @return Invisibly returns `relative_profile`.
#' @family hermosillo_analysis
#' @export
plot_profile_likelihood <- function(values,
                                    relative_profile,
                                    x_label,
                                    interval = c(lower = NA_real_, upper = NA_real_),
                                    interval_label = NULL,
                                    main = "Relative Profile Likelihood") {
  if (is.null(interval_label)) {
    interval_label <- format_interval_labels(interval)
  }

  graphics::plot(
    values,
    relative_profile,
    type = "l",
    lty = 1,
    lwd = 2,
    col = 1,
    cex.lab = 1.25,
    cex.axis = 1.25,
    main = main,
    xlab = x_label,
    ylim = c(0, 1),
    ylab = "Relative Profile Likelihood",
    cex.main = 1.30
  )

  if (!anyNA(interval)) {
    graphics::points(interval[["lower"]], 0.015, pch = -9658, col = 1, cex = 1)
    graphics::points(interval[["upper"]], 0.015, pch = -9668, col = 1, cex = 1)
  }

  if (!is.null(interval_label)) {
    graphics::legend(
      "topright",
      interval_label,
      pch = c(-9658, -9668)[seq_along(interval_label)],
      col = rep(1, length(interval_label)),
      bty = "n",
      title.adj = 0.35,
      y.intersp = 0.5,
      inset = -0.025,
      title = "95%CI"
    )
  }

  invisible(relative_profile)
}
