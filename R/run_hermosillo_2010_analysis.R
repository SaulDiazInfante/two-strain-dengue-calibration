#' Run the Hermosillo 2010 Analysis Workflow
#'
#' @description
#' Runs the full calibration and visualization workflow for the Hermosillo 2010
#' dengue case study using a compact nested configuration. The workflow reads
#' defaults from the packaged JSON config, supports JSON override files, and
#' also accepts in-memory config overrides for quick experiments.
#'
#' @param config Optional named list with nested sections such as `data`,
#'   `observed_rows`, `grids`, `outputs`, and `optimisation`.
#' @param config_path Optional path to a JSON file containing configuration
#'   overrides merged on top of the packaged defaults.
#' @param data Optional study data overriding the data source specified in the
#'   config.
#' @param ... Optional legacy flat arguments from the previous interface. These
#'   are converted into config overrides for backwards compatibility.
#'
#' @return A named list containing the normalized config, fitted parameters,
#'   likelihood surfaces, profile summaries, reproduction-number estimates,
#'   summary data frames, and any generated output file paths.
#'
#' @examples
#' analysis <- run_hermosillo_2010_analysis(
#'     config = list(
#'         observed_rows = 1:3,
#'         grids = list(
#'             time = list(
#'                 full = seq(0, 5, by = 0.5),
#'                 calibration = seq(0, 2, by = 0.5),
#'                 fit = seq(0, 5, by = 0.5)
#'             ),
#'             parameters = list(
#'                 beta_h = seq(1, 1.4, length.out = 3),
#'                 beta_m = seq(0.001, 0.003, length.out = 3),
#'                 r01 = seq(1.7, 1.9, length.out = 3),
#'                 r02 = seq(0.003, 0.0034, length.out = 3),
#'                 r0 = seq(1.3, 1.4, length.out = 3)
#'             )
#'         ),
#'         outputs = list(write_outputs = FALSE, show_progress = FALSE)
#'     )
#' )
#'
#' analysis$reproduction_numbers
#' @family hermosillo_analysis
#' @export
run_hermosillo_2010_analysis <- function(config = NULL,
                                         config_path = NULL,
                                         data = NULL,
                                         ...) {
  context <- prepare_hermosillo_2010_analysis_context(
    config = config,
    config_path = config_path,
    data = data,
    legacy_arguments = list(...)
  )
  progress <- create_analysis_progress(
    enabled = context$outputs$show_progress,
    total_steps = 5L +
      nrow(context$surface_table) +
      as.integer(context$outputs$write_outputs)
  )

  progress$tick("Running baseline simulation")
  baseline_output <- run_model(
    context$initial_state,
    context$time_grids$full,
    context$all_parameters
  )

  primary_surface <- evaluate_primary_hermosillo_2010_surface(context, progress)
  mle <- optimise_hermosillo_2010_parameters(context, progress)
  primary_profiles <- build_hermosillo_2010_primary_profiles(
    parameter_grids = context$parameter_grids,
    global_loglikelihood = primary_surface$loglikelihood
  )
  reproduction_numbers <- compute_hermosillo_2010_reproduction_numbers(
    mle_parameters = mle$parameters,
    initial_state = context$initial_state,
    fixed_parameters = context$fixed_parameters
  )
  reparameterised_surface_results <-
    evaluate_reparameterised_hermosillo_2010_surfaces(context, progress)
  reparameterised_profiles <- build_hermosillo_2010_reparameterised_profiles(
    surface_table = context$surface_table,
    parameter_grids = context$parameter_grids,
    reparameterised_surfaces = reparameterised_surface_results$surfaces
  )
  profile_summary <- build_hermosillo_2010_profile_summary(
    surface_table = context$surface_table,
    primary_profiles = primary_profiles,
    reparameterised_profiles = reparameterised_profiles
  )

  output_files <- list(
    figures = list(),
    calibration = c(
      primary_surface$calibration_files,
      reparameterised_surface_results$calibration_files
    )
  )

  if (context$outputs$write_outputs) {
    progress$tick("Writing figures and calibration tables")
    output_files$figures <- write_hermosillo_2010_figures(
      context = context,
      baseline_output = baseline_output,
      primary_surface = primary_surface,
      mle = mle,
      primary_profiles = primary_profiles,
      reproduction_numbers = reproduction_numbers,
      reparameterised_surfaces = reparameterised_surface_results$surfaces,
      reparameterised_profiles = reparameterised_profiles
    )
  }

  progress$tick("Compiling analysis results")
  build_hermosillo_2010_analysis_result(
    context = context,
    baseline_output = baseline_output,
    mle = mle,
    reproduction_numbers = reproduction_numbers,
    primary_surface = primary_surface,
    primary_profiles = primary_profiles,
    reparameterised_surfaces = reparameterised_surface_results$surfaces,
    reparameterised_profiles = reparameterised_profiles,
    profile_summary = profile_summary,
    output_files = output_files
  )
}
