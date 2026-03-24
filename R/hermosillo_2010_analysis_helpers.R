#' Prepare the Hermosillo 2010 Analysis Context
#'
#' @param config Optional nested config overrides.
#' @param config_path Optional path to a JSON config override file.
#' @param data Optional data frame overriding the configured data source.
#' @param legacy_arguments Named list of legacy flat arguments.
#'
#' @return A list of normalized workflow inputs derived from the config.
#' @keywords internal
#' @noRd
prepare_hermosillo_2010_analysis_context <- function(config,
                                                     config_path,
                                                     data,
                                                     legacy_arguments) {
  legacy_config <- if (length(legacy_arguments) > 0L) {
    legacy_hermosillo_2010_analysis_overrides(legacy_arguments)
  } else {
    list()
  }

  analysis_config <- load_hermosillo_2010_analysis_config(path = config_path)
  analysis_config <- normalise_hermosillo_2010_analysis_config(
    merge_nested_lists(
      analysis_config,
      merge_nested_lists(config, legacy_config)
    )
  )

  if (is.null(data)) {
    data_path <- analysis_config$data$path
    data <- if (is.null(data_path)) {
      load_study_data()
    } else {
      load_study_data(path = data_path)
    }
  }

  observed_rows <- analysis_config$observed_rows
  time_grids <- analysis_config$grids$time
  parameter_grids <- analysis_config$grids$parameters
  outputs <- analysis_config$outputs
  optimisation <- analysis_config$optimisation
  surface_table <- build_reparameterised_surface_table()
  timestamp <- outputs$timestamp

  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y_%m_%d__%H_%M")
  }

  all_parameters <- default_parameters()
  fixed_parameters <- all_parameters[fixed_parameter_names]
  initial_state <- default_state(all_parameters)
  full_rows <- seq_len(nrow(data))
  week_origin <- min(data$week[observed_rows])

  if (is.null(time_grids$fit)) {
    time_grids$fit <- seq(0, max(data$week) - week_origin, by = 0.1)
  }

  observation_times <- data$week[observed_rows] - week_origin
  observed_cases <- as.matrix(
    data[
      observed_rows,
      c("classical", "hemorrhagic")
    ]
  )

  list(
    analysis_config = analysis_config,
    timestamp = timestamp,
    data = data,
    observed_rows = observed_rows,
    full_rows = full_rows,
    week_origin = week_origin,
    observation_times = observation_times,
    observed_cases = observed_cases,
    time_grids = time_grids,
    parameter_grids = parameter_grids,
    outputs = outputs,
    optimisation = optimisation,
    surface_table = surface_table,
    all_parameters = all_parameters,
    fixed_parameters = fixed_parameters,
    initial_state = initial_state
  )
}

#' Evaluate the Primary Hermosillo Likelihood Surface
#'
#' @param context Analysis context built by
#'   [prepare_hermosillo_2010_analysis_context()].
#' @param progress Stage reporter built by [create_analysis_progress()].
#'
#' @return A list with the global likelihood surface, its relative version, and
#'   any written calibration file paths.
#' @keywords internal
#' @noRd
evaluate_primary_hermosillo_2010_surface <- function(context, progress) {
  parameter_grids <- context$parameter_grids

  progress$tick("Evaluating betaH / betaM likelihood surface")
  global_loglikelihood <- evaluate_loglikelihood_grid_with_progress(
    parameter_grids$beta_h,
    parameter_grids$beta_m,
    function(par) {
      negative_loglikelihood_original_scale(
        par = par,
        observation_times = context$observation_times,
        observed_cases = context$observed_cases,
        initial_state = context$initial_state,
        time_grid = context$time_grids$calibration,
        fixed_parameters = context$fixed_parameters
      )
    },
    show_progress = context$outputs$show_progress,
    progress_label = sprintf(
      "  betaH / betaM grid (%d rows x %d columns)",
      length(parameter_grids$beta_h),
      length(parameter_grids$beta_m)
    )
  )

  calibration_files <- list()

  if (context$outputs$write_outputs) {
    calibration_files$beta_h_beta_m <- build_hermosillo_2010_calibration_files(
      outputs = context$outputs,
      timestamp = context$timestamp,
      log_stem = "loglik_surface_beta_h_beta_m_poisson",
      relative_stem = "relative_loglik_surface_beta_h_beta_m_poisson"
    )
    relative_loglikelihood <- write_likelihood_outputs(
      global_loglikelihood,
      calibration_files$beta_h_beta_m$loglikelihood,
      calibration_files$beta_h_beta_m$relative
    )
  } else {
    relative_loglikelihood <- relative_likelihood(global_loglikelihood)
  }

  list(
    loglikelihood = global_loglikelihood,
    relative = relative_loglikelihood,
    calibration_files = calibration_files
  )
}

#' Optimise Hermosillo Transmission Parameters
#'
#' @param context Analysis context built by
#'   [prepare_hermosillo_2010_analysis_context()].
#' @param progress Stage reporter built by [create_analysis_progress()].
#'
#' @return A list containing the optimiser fit, fitted parameters, and the
#'   model trajectory at the fitted values.
#' @keywords internal
#' @noRd
optimise_hermosillo_2010_parameters <- function(context, progress) {
  progress$tick("Optimising betaH and betaM")
  mle_fit <- stats::constrOptim(
    theta = context$optimisation$start,
    f = negative_loglikelihood_original_scale,
    grad = NULL,
    ui = diag(2),
    ci = c(0, 0),
    mu = 1e-04,
    method = "Nelder-Mead",
    outer.iterations = 100,
    outer.eps = 1e-05,
    observation_times = context$observation_times,
    observed_cases = context$observed_cases,
    initial_state = context$initial_state,
    time_grid = context$time_grids$calibration,
    fixed_parameters = context$fixed_parameters
  )
  mle_parameters <- stats::setNames(mle_fit$par, c("betaH", "betaM"))

  progress$tick("Running fitted model trajectory")
  fit_output <- run_model(
    context$initial_state,
    context$time_grids$fit,
    build_parameter_vector(
      context$fixed_parameters,
      mle_parameters[["betaH"]],
      mle_parameters[["betaM"]]
    )
  )

  list(
    optimisation = mle_fit,
    parameters = mle_parameters,
    fit_output = fit_output
  )
}

#' Build Primary Hermosillo Profile Likelihoods
#'
#' @param parameter_grids Nested parameter grids from the normalized config.
#' @param global_loglikelihood Matrix of `betaH` by `betaM` log-likelihoods.
#'
#' @return A named list containing the primary profile likelihoods and
#'   intervals for `beta_h` and `beta_m`.
#' @keywords internal
#' @noRd
build_hermosillo_2010_primary_profiles <- function(parameter_grids,
                                                   global_loglikelihood) {
  beta_h_relative <- relative_likelihood(apply(global_loglikelihood, 1, max))
  beta_m_relative <- relative_likelihood(apply(global_loglikelihood, 2, max))

  list(
    beta_h = list(
      values = parameter_grids$beta_h,
      relative = beta_h_relative,
      interval = profile_interval(parameter_grids$beta_h, beta_h_relative)
    ),
    beta_m = list(
      values = parameter_grids$beta_m,
      relative = beta_m_relative,
      interval = profile_interval(parameter_grids$beta_m, beta_m_relative)
    )
  )
}

#' Compute Hermosillo Reproduction Numbers
#'
#' @param mle_parameters Named numeric vector containing `betaH` and `betaM`.
#' @param initial_state Named numeric state vector used by [run_model()].
#' @param fixed_parameters Named numeric vector of fixed model parameters.
#'
#' @return A named numeric vector containing `R01`, `R02`, and `R0`.
#' @keywords internal
#' @noRd
compute_hermosillo_2010_reproduction_numbers <- function(mle_parameters,
                                                         initial_state,
                                                         fixed_parameters) {
  reproduction_components <- compute_reproduction_components(
    initial_state,
    fixed_parameters
  )
  pi_r <- mle_parameters[["betaH"]] *
    mle_parameters[["betaM"]] *
    reproduction_components[["C1"]]

  c(
    R01 = pi_r * reproduction_components[["C2"]],
    R02 = pi_r * reproduction_components[["C3"]],
    R0 = sqrt((pi_r * reproduction_components[["C2"]]) +
      (pi_r * reproduction_components[["C3"]]))
  )
}

#' Evaluate Reparameterised Hermosillo Likelihood Surfaces
#'
#' @param context Analysis context built by
#'   [prepare_hermosillo_2010_analysis_context()].
#' @param progress Stage reporter built by [create_analysis_progress()].
#'
#' @return A list containing the reparameterised surfaces and any calibration
#'   files written while evaluating them.
#' @keywords internal
#' @noRd
evaluate_reparameterised_hermosillo_2010_surfaces <- function(context,
                                                              progress) {
  reparameterised_surfaces <- setNames(
    vector("list", nrow(context$surface_table)),
    context$surface_table$surface_id
  )
  calibration_files <- list()

  for (i in seq_len(nrow(context$surface_table))) {
    spec <- context$surface_table[i, , drop = FALSE]
    surface_id <- spec$surface_id[[1]]
    grid_values <- context$parameter_grids[[spec$grid_key[[1]]]]
    surface_files <- if (context$outputs$write_outputs) {
      build_hermosillo_2010_calibration_files(
        outputs = context$outputs,
        timestamp = context$timestamp,
        log_stem = spec$calibration_log_stem[[1]],
        relative_stem = spec$calibration_relative_stem[[1]]
      )
    } else {
      list(loglikelihood = NULL, relative = NULL)
    }

    progress$tick(sprintf(
      "Evaluating betaH / %s likelihood surface",
      spec$display_label[[1]]
    ))
    reparameterised_surfaces[[surface_id]] <- evaluate_surface_for_mode(
      beta_values = context$parameter_grids$beta_h,
      reproduction_values = grid_values,
      mode = spec$mode[[1]],
      observation_times = context$observation_times,
      observed_cases = context$observed_cases,
      initial_state = context$initial_state,
      time_grid = context$time_grids$calibration,
      fixed_parameters = context$fixed_parameters,
      show_progress = context$outputs$show_progress,
      progress_label = sprintf(
        "  betaH / %s grid (%d rows x %d columns)",
        spec$display_label[[1]],
        length(context$parameter_grids$beta_h),
        length(grid_values)
      ),
      log_file = surface_files$loglikelihood,
      relative_file = surface_files$relative
    )

    if (context$outputs$write_outputs) {
      calibration_files[[surface_id]] <- reparameterised_surfaces[[surface_id]]$files
    }
  }

  list(
    surfaces = reparameterised_surfaces,
    calibration_files = calibration_files
  )
}

#' Build Reparameterised Hermosillo Profiles
#'
#' @param surface_table Metadata table built by
#'   [build_reparameterised_surface_table()].
#' @param parameter_grids Nested parameter grids from the normalized config.
#' @param reparameterised_surfaces Named list of likelihood surfaces.
#'
#' @return A named list of profile likelihoods and intervals for each
#'   reparameterised surface.
#' @keywords internal
#' @noRd
build_hermosillo_2010_reparameterised_profiles <- function(surface_table,
                                                           parameter_grids,
                                                           reparameterised_surfaces) {
  reparameterised_profiles <- setNames(
    vector("list", nrow(surface_table)),
    surface_table$surface_id
  )

  for (i in seq_len(nrow(surface_table))) {
    spec <- surface_table[i, , drop = FALSE]
    surface_id <- spec$surface_id[[1]]
    grid_values <- parameter_grids[[spec$grid_key[[1]]]]
    relative_profile <- relative_likelihood(
      apply(reparameterised_surfaces[[surface_id]]$loglikelihood, 2, max)
    )

    reparameterised_profiles[[surface_id]] <- list(
      values = grid_values,
      relative = relative_profile,
      interval = profile_interval(grid_values, relative_profile)
    )
  }

  reparameterised_profiles
}

#' Build the Hermosillo Profile Summary Table
#'
#' @param surface_table Metadata table built by
#'   [build_reparameterised_surface_table()].
#' @param primary_profiles Named list of `beta_h` and `beta_m` profiles.
#' @param reparameterised_profiles Named list of reproduction-number profiles.
#'
#' @return A data frame with lower and upper bounds for each profile interval.
#' @keywords internal
#' @noRd
build_hermosillo_2010_profile_summary <- function(surface_table,
                                                  primary_profiles,
                                                  reparameterised_profiles) {
  data.frame(
    parameter = c("beta_h", "beta_m", surface_table$surface_id),
    lower = c(
      primary_profiles$beta_h$interval[["lower"]],
      primary_profiles$beta_m$interval[["lower"]],
      vapply(
        reparameterised_profiles,
        function(profile) profile$interval[["lower"]],
        numeric(1)
      )
    ),
    upper = c(
      primary_profiles$beta_h$interval[["upper"]],
      primary_profiles$beta_m$interval[["upper"]],
      vapply(
        reparameterised_profiles,
        function(profile) profile$interval[["upper"]],
        numeric(1)
      )
    ),
    stringsAsFactors = FALSE
  )
}

#' Write Hermosillo Analysis Figures
#'
#' @param context Analysis context built by
#'   [prepare_hermosillo_2010_analysis_context()].
#' @param baseline_output Model trajectory at the default parameters.
#' @param primary_surface Primary `betaH` by `betaM` likelihood surface.
#' @param mle Fitted-parameter results returned by
#'   [optimise_hermosillo_2010_parameters()].
#' @param primary_profiles Named list of `beta_h` and `beta_m` profiles.
#' @param reproduction_numbers Named numeric vector containing `R01`, `R02`,
#'   and `R0`.
#' @param reparameterised_surfaces Named list of reproduction-number surfaces.
#' @param reparameterised_profiles Named list of reproduction-number profiles.
#'
#' @return A named list of written figure paths.
#' @keywords internal
#' @noRd
write_hermosillo_2010_figures <- function(context,
                                          baseline_output,
                                          primary_surface,
                                          mle,
                                          primary_profiles,
                                          reproduction_numbers,
                                          reparameterised_surfaces,
                                          reparameterised_profiles) {
  c(
    write_primary_hermosillo_2010_figures(
      context = context,
      baseline_output = baseline_output,
      primary_surface = primary_surface,
      mle = mle,
      primary_profiles = primary_profiles
    ),
    write_reparameterised_hermosillo_2010_figures(
      context = context,
      reparameterised_surfaces = reparameterised_surfaces,
      reparameterised_profiles = reparameterised_profiles,
      mle_parameters = mle$parameters,
      reproduction_numbers = reproduction_numbers
    )
  )
}

#' Build the Hermosillo Analysis Result Object
#'
#' @param context Analysis context built by
#'   [prepare_hermosillo_2010_analysis_context()].
#' @param baseline_output Model trajectory at the default parameters.
#' @param mle Fitted-parameter results returned by
#'   [optimise_hermosillo_2010_parameters()].
#' @param reproduction_numbers Named numeric vector containing `R01`, `R02`,
#'   and `R0`.
#' @param primary_surface Primary `betaH` by `betaM` likelihood surface.
#' @param primary_profiles Named list of `beta_h` and `beta_m` profiles.
#' @param reparameterised_surfaces Named list of reproduction-number surfaces.
#' @param reparameterised_profiles Named list of reproduction-number profiles.
#' @param profile_summary Data frame summarizing profile intervals.
#' @param output_files Named list of generated calibration and figure files.
#'
#' @return The structured result returned by
#'   [run_hermosillo_2010_analysis()].
#' @keywords internal
#' @noRd
build_hermosillo_2010_analysis_result <- function(context,
                                                  baseline_output,
                                                  mle,
                                                  reproduction_numbers,
                                                  primary_surface,
                                                  primary_profiles,
                                                  reparameterised_surfaces,
                                                  reparameterised_profiles,
                                                  profile_summary,
                                                  output_files) {
  list(
    config = context$analysis_config,
    timestamp = context$timestamp,
    week_origin = context$week_origin,
    data = context$data,
    observed_rows = context$observed_rows,
    initial_state = context$initial_state,
    fixed_parameters = context$fixed_parameters,
    baseline_output = baseline_output,
    mle = mle,
    reproduction_numbers = reproduction_numbers,
    surface_metadata = context$surface_table[
      c(
        "surface_id",
        "mode",
        "grid_key",
        "number_key",
        "display_label",
        "calibration_log_stem",
        "calibration_relative_stem",
        "surface_figure_stem",
        "profile_figure_stem"
      )
    ],
    surfaces = c(
      list(
        beta_h_beta_m = list(
          loglikelihood = primary_surface$loglikelihood,
          relative = primary_surface$relative
        )
      ),
      reparameterised_surfaces
    ),
    profiles = c(primary_profiles, reparameterised_profiles),
    profile_summary = profile_summary,
    output_files = output_files
  )
}

build_hermosillo_2010_calibration_files <- function(outputs,
                                                    timestamp,
                                                    log_stem,
                                                    relative_stem) {
  list(
    loglikelihood = build_output_path(
      outputs$dir,
      "calibration",
      log_stem,
      "csv",
      timestamp
    ),
    relative = build_output_path(
      outputs$dir,
      "calibration",
      relative_stem,
      "csv",
      timestamp
    )
  )
}

write_primary_hermosillo_2010_figures <- function(context,
                                                  baseline_output,
                                                  primary_surface,
                                                  mle,
                                                  primary_profiles) {
  data <- context$data
  full_rows <- context$full_rows
  observed_rows <- context$observed_rows
  outputs <- context$outputs
  parameter_grids <- context$parameter_grids
  timestamp <- context$timestamp
  week_origin <- context$week_origin

  list(
    observed_cases_full = save_plot(
      build_output_path(outputs$dir, "figures", "observed_cases_full", "jpeg", timestamp),
      function() {
        plot_observed_cases(
          data = data,
          rows = full_rows,
          title = "Observed Cases: Hermosillo 2010",
          ylim = c(0, 200)
        )
      }
    ),
    model_fit_full_window = save_plot(
      build_output_path(outputs$dir, "figures", "model_fit_full_window", "jpeg", timestamp),
      function() {
        plot_model_fit(
          model_output = baseline_output,
          data = data,
          rows = full_rows,
          title = "Baseline Model Fit: Hermosillo 2010",
          week_origin = week_origin
        )
      }
    ),
    observed_cases_calibration_window = save_plot(
      build_output_path(
        outputs$dir,
        "figures",
        "observed_cases_calibration_window",
        "jpeg",
        timestamp
      ),
      function() {
        plot_observed_cases(
          data = data,
          rows = observed_rows,
          title = "Observed Cases: Calibration Window",
          ylim = c(0, 150)
        )
      }
    ),
    likelihood_surface_beta_h_beta_m = save_plot(
      build_output_path(
        outputs$dir,
        "figures",
        "likelihood_surface_beta_h_beta_m",
        "jpeg",
        timestamp
      ),
      function() {
        plot_likelihood_surface(
          parameter_grids$beta_h,
          parameter_grids$beta_m,
          primary_surface$relative,
          x_label = expression(beta[H]),
          y_label = expression(beta[M]),
          point = unname(mle$parameters),
          label = sprintf(
            "(betaH, betaM) = (%.4f, %.4f)",
            mle$parameters[["betaH"]],
            mle$parameters[["betaM"]]
          )
        )
      }
    ),
    model_fit_calibration_window = save_plot(
      build_output_path(
        outputs$dir,
        "figures",
        "model_fit_calibration_window",
        "jpeg",
        timestamp
      ),
      function() {
        plot_model_fit(
          model_output = mle$fit_output,
          data = data,
          rows = observed_rows,
          title = "Poisson-Strain Model Fit: Calibration Window",
          week_origin = week_origin,
          xlim = range(data$week[observed_rows])
        )
      }
    ),
    model_fit_full_data = save_plot(
      build_output_path(outputs$dir, "figures", "model_fit_full_data", "jpeg", timestamp),
      function() {
        plot_model_fit(
          model_output = mle$fit_output,
          data = data,
          rows = full_rows,
          title = "Poisson-Strain Model Fit: Full Study Window",
          week_origin = week_origin,
          xlim = range(data$week[full_rows])
        )
      }
    ),
    profile_likelihood_beta_h = save_plot(
      build_output_path(outputs$dir, "figures", "profile_likelihood_beta_h", "jpeg", timestamp),
      function() {
        plot_profile_likelihood(
          parameter_grids$beta_h,
          primary_profiles$beta_h$relative,
          x_label = expression(beta[H]),
          interval = primary_profiles$beta_h$interval,
          main = expression(paste("Profile Likelihood for ", beta[H]))
        )
      }
    ),
    profile_likelihood_beta_m = save_plot(
      build_output_path(outputs$dir, "figures", "profile_likelihood_beta_m", "jpeg", timestamp),
      function() {
        plot_profile_likelihood(
          parameter_grids$beta_m,
          primary_profiles$beta_m$relative,
          x_label = expression(beta[M]),
          interval = primary_profiles$beta_m$interval,
          main = expression(paste("Profile Likelihood for ", beta[M]))
        )
      }
    )
  )
}

write_reparameterised_hermosillo_2010_figures <- function(context,
                                                          reparameterised_surfaces,
                                                          reparameterised_profiles,
                                                          mle_parameters,
                                                          reproduction_numbers) {
  figure_files <- list()

  for (i in seq_len(nrow(context$surface_table))) {
    spec <- context$surface_table[i, , drop = FALSE]
    surface_id <- spec$surface_id[[1]]
    grid_values <- context$parameter_grids[[spec$grid_key[[1]]]]
    profile <- reparameterised_profiles[[surface_id]]

    figure_files[[spec$surface_figure_stem[[1]]]] <- save_plot(
      build_output_path(
        context$outputs$dir,
        "figures",
        spec$surface_figure_stem[[1]],
        "jpeg",
        context$timestamp
      ),
      function() {
        plot_likelihood_surface(
          context$parameter_grids$beta_h,
          grid_values,
          reparameterised_surfaces[[surface_id]]$relative,
          x_label = expression(beta[H]),
          y_label = spec$y_axis_label[[1]],
          point = c(
            mle_parameters[["betaH"]],
            reproduction_numbers[[spec$number_key[[1]]]]
          ),
          label = sprintf(
            "(betaH, %s) = (%.4f, %.4f)",
            spec$display_label[[1]],
            mle_parameters[["betaH"]],
            reproduction_numbers[[spec$number_key[[1]]]]
          )
        )
      }
    )

    figure_files[[spec$profile_figure_stem[[1]]]] <- save_plot(
      build_output_path(
        context$outputs$dir,
        "figures",
        spec$profile_figure_stem[[1]],
        "jpeg",
        context$timestamp
      ),
      function() {
        plot_profile_likelihood(
          grid_values,
          profile$relative,
          x_label = spec$y_axis_label[[1]],
          interval = profile$interval,
          main = spec$profile_title[[1]]
        )
      }
    )
  }

  figure_files
}
