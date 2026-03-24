#' Build the Reparameterised Surface Metadata Table
#'
#' @return A data frame describing the reproduction-number likelihood surfaces.
#' @keywords internal
#' @noRd
build_reparameterised_surface_table <- function() {
  data.frame(
    surface_id = c("r01", "r02", "r0"),
    mode = c("R01", "R02", "R0"),
    grid_key = c("r01", "r02", "r0"),
    number_key = c("R01", "R02", "R0"),
    display_label = c("R1", "R2", "R0"),
    calibration_log_stem = c(
      "loglik_surface_beta_h_r1_poisson",
      "loglik_surface_beta_h_r2_poisson",
      "loglik_surface_beta_h_r0_poisson"
    ),
    calibration_relative_stem = c(
      "relative_loglik_surface_beta_h_r1_poisson",
      "relative_loglik_surface_beta_h_r2_poisson",
      "relative_loglik_surface_beta_h_r0_poisson"
    ),
    surface_figure_stem = c(
      "likelihood_surface_beta_h_r1",
      "likelihood_surface_beta_h_r2",
      "likelihood_surface_beta_h_r0"
    ),
    profile_figure_stem = c(
      "profile_likelihood_r1",
      "profile_likelihood_r2",
      "profile_likelihood_r0"
    ),
    y_axis_label = I(list(
      expression(R[1]),
      expression(R[2]),
      expression(R[0])
    )),
    profile_title = I(list(
      expression(paste("Profile Likelihood for ", R[1])),
      expression(paste("Profile Likelihood for ", R[2])),
      expression(paste("Profile Likelihood for ", R[0]))
    )),
    stringsAsFactors = FALSE
  )
}
