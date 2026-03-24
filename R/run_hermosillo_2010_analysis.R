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
    legacy_arguments <- list(...)
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
    progress <- create_analysis_progress(
        enabled = outputs$show_progress,
        total_steps = 5L + nrow(surface_table) + as.integer(outputs$write_outputs)
    )

    progress$tick("Running baseline simulation")
    baseline_output <- run_model(initial_state, time_grids$full, all_parameters)

    progress$tick("Evaluating betaH / betaM likelihood surface")
    global_loglikelihood <- evaluate_loglikelihood_grid_with_progress(
        parameter_grids$beta_h,
        parameter_grids$beta_m,
        function(par) {
            negative_loglikelihood_original_scale(
                par = par,
                observation_times = observation_times,
                observed_cases = observed_cases,
                initial_state = initial_state,
                time_grid = time_grids$calibration,
                fixed_parameters = fixed_parameters
            )
        },
        show_progress = outputs$show_progress,
        progress_label = sprintf(
            "  betaH / betaM grid (%d rows x %d columns)",
            length(parameter_grids$beta_h),
            length(parameter_grids$beta_m)
        )
    )

    output_files <- list(figures = list(), calibration = list())

    if (outputs$write_outputs) {
        beta_h_beta_m_log <- build_output_path(
            outputs$dir,
            "calibration",
            "loglik_surface_beta_h_beta_m_poisson",
            "csv",
            timestamp
        )
        beta_h_beta_m_relative <- build_output_path(
            outputs$dir,
            "calibration",
            "relative_loglik_surface_beta_h_beta_m_poisson",
            "csv",
            timestamp
        )
        relative_global_loglikelihood <- write_likelihood_outputs(
            global_loglikelihood,
            beta_h_beta_m_log,
            beta_h_beta_m_relative
        )
        output_files$calibration$beta_h_beta_m <- list(
            loglikelihood = beta_h_beta_m_log,
            relative = beta_h_beta_m_relative
        )
    } else {
        relative_global_loglikelihood <- relative_likelihood(global_loglikelihood)
    }

    progress$tick("Optimising betaH and betaM")
    mle_fit <- stats::constrOptim(
        theta = optimisation$start,
        f = negative_loglikelihood_original_scale,
        grad = NULL,
        ui = diag(2),
        ci = c(0, 0),
        mu = 1e-04,
        method = "Nelder-Mead",
        outer.iterations = 100,
        outer.eps = 1e-05,
        observation_times = observation_times,
        observed_cases = observed_cases,
        initial_state = initial_state,
        time_grid = time_grids$calibration,
        fixed_parameters = fixed_parameters
    )
    mle_parameters <- stats::setNames(mle_fit$par, c("betaH", "betaM"))

    progress$tick("Running fitted model trajectory")
    fit_output <- run_model(
        initial_state,
        time_grids$fit,
        build_parameter_vector(
            fixed_parameters,
            mle_parameters[["betaH"]],
            mle_parameters[["betaM"]]
        )
    )

    beta_h_profile <- relative_likelihood(apply(global_loglikelihood, 1, max))
    beta_h_interval <- profile_interval(parameter_grids$beta_h, beta_h_profile)
    beta_m_profile <- relative_likelihood(apply(global_loglikelihood, 2, max))
    beta_m_interval <- profile_interval(parameter_grids$beta_m, beta_m_profile)

    reproduction_components <- compute_reproduction_components(
        initial_state,
        fixed_parameters
    )
    pi_r <- mle_parameters[["betaH"]] *
        mle_parameters[["betaM"]] *
        reproduction_components[["C1"]]
    reproduction_numbers <- c(
        R01 = pi_r * reproduction_components[["C2"]],
        R02 = pi_r * reproduction_components[["C3"]],
        R0 = sqrt((pi_r * reproduction_components[["C2"]]) +
            (pi_r * reproduction_components[["C3"]]))
    )

    reparameterised_surfaces <- setNames(
        vector("list", nrow(surface_table)),
        surface_table$surface_id
    )

    for (i in seq_len(nrow(surface_table))) {
        spec <- surface_table[i, , drop = FALSE]
        grid_values <- parameter_grids[[spec$grid_key[[1]]]]

        progress$tick(sprintf(
            "Evaluating betaH / %s likelihood surface",
            spec$display_label[[1]]
        ))
        reparameterised_surfaces[[spec$surface_id[[1]]]] <- evaluate_surface_for_mode(
            beta_values = parameter_grids$beta_h,
            reproduction_values = grid_values,
            mode = spec$mode[[1]],
            observation_times = observation_times,
            observed_cases = observed_cases,
            initial_state = initial_state,
            time_grid = time_grids$calibration,
            fixed_parameters = fixed_parameters,
            show_progress = outputs$show_progress,
            progress_label = sprintf(
                "  betaH / %s grid (%d rows x %d columns)",
                spec$display_label[[1]],
                length(parameter_grids$beta_h),
                length(grid_values)
            ),
            log_file = if (outputs$write_outputs) {
                build_output_path(
                    outputs$dir,
                    "calibration",
                    spec$calibration_log_stem[[1]],
                    "csv",
                    timestamp
                )
            } else {
                NULL
            },
            relative_file = if (outputs$write_outputs) {
                build_output_path(
                    outputs$dir,
                    "calibration",
                    spec$calibration_relative_stem[[1]],
                    "csv",
                    timestamp
                )
            } else {
                NULL
            }
        )
    }

    if (outputs$write_outputs) {
        for (surface_id in names(reparameterised_surfaces)) {
            output_files$calibration[[surface_id]] <- reparameterised_surfaces[[surface_id]]$files
        }
    }

    reparameterised_profiles <- setNames(
        vector("list", nrow(surface_table)),
        surface_table$surface_id
    )

    for (i in seq_len(nrow(surface_table))) {
        spec <- surface_table[i, , drop = FALSE]
        grid_values <- parameter_grids[[spec$grid_key[[1]]]]
        relative_profile <- relative_likelihood(
            apply(reparameterised_surfaces[[spec$surface_id[[1]]]]$loglikelihood, 2, max)
        )

        reparameterised_profiles[[spec$surface_id[[1]]]] <- list(
            values = grid_values,
            relative = relative_profile,
            interval = profile_interval(grid_values, relative_profile)
        )
    }

    profile_summary <- data.frame(
        parameter = c("beta_h", "beta_m", surface_table$surface_id),
        lower = c(
            beta_h_interval[["lower"]],
            beta_m_interval[["lower"]],
            vapply(
                reparameterised_profiles,
                function(profile) profile$interval[["lower"]],
                numeric(1)
            )
        ),
        upper = c(
            beta_h_interval[["upper"]],
            beta_m_interval[["upper"]],
            vapply(
                reparameterised_profiles,
                function(profile) profile$interval[["upper"]],
                numeric(1)
            )
        ),
        stringsAsFactors = FALSE
    )

    if (outputs$write_outputs) {
        progress$tick("Writing figures and calibration tables")
        output_files$figures$observed_cases_full <- save_plot(
            build_output_path(outputs$dir, "figures", "observed_cases_full", "jpeg", timestamp),
            function() {
                plot_observed_cases(
                    data = data,
                    rows = full_rows,
                    title = "Observed Cases: Hermosillo 2010",
                    ylim = c(0, 200)
                )
            }
        )
        output_files$figures$model_fit_full_window <- save_plot(
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
        )
        output_files$figures$observed_cases_calibration_window <- save_plot(
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
        )
        output_files$figures$likelihood_surface_beta_h_beta_m <- save_plot(
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
                    relative_global_loglikelihood,
                    x_label = expression(beta[H]),
                    y_label = expression(beta[M]),
                    point = unname(mle_parameters),
                    label = sprintf(
                        "(betaH, betaM) = (%.4f, %.4f)",
                        mle_parameters[["betaH"]],
                        mle_parameters[["betaM"]]
                    )
                )
            }
        )
        output_files$figures$model_fit_calibration_window <- save_plot(
            build_output_path(
                outputs$dir,
                "figures",
                "model_fit_calibration_window",
                "jpeg",
                timestamp
            ),
            function() {
                plot_model_fit(
                    model_output = fit_output,
                    data = data,
                    rows = observed_rows,
                    title = "Poisson-Strain Model Fit: Calibration Window",
                    week_origin = week_origin,
                    xlim = range(data$week[observed_rows])
                )
            }
        )
        output_files$figures$model_fit_full_data <- save_plot(
            build_output_path(outputs$dir, "figures", "model_fit_full_data", "jpeg", timestamp),
            function() {
                plot_model_fit(
                    model_output = fit_output,
                    data = data,
                    rows = full_rows,
                    title = "Poisson-Strain Model Fit: Full Study Window",
                    week_origin = week_origin,
                    xlim = range(data$week[full_rows])
                )
            }
        )
        output_files$figures$profile_likelihood_beta_h <- save_plot(
            build_output_path(outputs$dir, "figures", "profile_likelihood_beta_h", "jpeg", timestamp),
            function() {
                plot_profile_likelihood(
                    parameter_grids$beta_h,
                    beta_h_profile,
                    x_label = expression(beta[H]),
                    interval = beta_h_interval,
                    main = expression(paste("Profile Likelihood for ", beta[H]))
                )
            }
        )
        output_files$figures$profile_likelihood_beta_m <- save_plot(
            build_output_path(outputs$dir, "figures", "profile_likelihood_beta_m", "jpeg", timestamp),
            function() {
                plot_profile_likelihood(
                    parameter_grids$beta_m,
                    beta_m_profile,
                    x_label = expression(beta[M]),
                    interval = beta_m_interval,
                    main = expression(paste("Profile Likelihood for ", beta[M]))
                )
            }
        )

        for (i in seq_len(nrow(surface_table))) {
            spec <- surface_table[i, , drop = FALSE]
            surface_id <- spec$surface_id[[1]]
            grid_values <- parameter_grids[[spec$grid_key[[1]]]]
            profile <- reparameterised_profiles[[surface_id]]

            output_files$figures[[spec$surface_figure_stem[[1]]]] <- save_plot(
                build_output_path(
                    outputs$dir,
                    "figures",
                    spec$surface_figure_stem[[1]],
                    "jpeg",
                    timestamp
                ),
                function() {
                    plot_likelihood_surface(
                        parameter_grids$beta_h,
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

            output_files$figures[[spec$profile_figure_stem[[1]]]] <- save_plot(
                build_output_path(
                    outputs$dir,
                    "figures",
                    spec$profile_figure_stem[[1]],
                    "jpeg",
                    timestamp
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
    }

    progress$tick("Compiling analysis results")
    list(
        config = analysis_config,
        timestamp = timestamp,
        week_origin = week_origin,
        data = data,
        observed_rows = observed_rows,
        initial_state = initial_state,
        fixed_parameters = fixed_parameters,
        baseline_output = baseline_output,
        mle = list(
            optimisation = mle_fit,
            parameters = mle_parameters,
            fit_output = fit_output
        ),
        reproduction_numbers = reproduction_numbers,
        surface_metadata = surface_table[
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
                    loglikelihood = global_loglikelihood,
                    relative = relative_global_loglikelihood
                )
            ),
            reparameterised_surfaces
        ),
        profiles = c(
            list(
                beta_h = list(
                    values = parameter_grids$beta_h,
                    relative = beta_h_profile,
                    interval = beta_h_interval
                ),
                beta_m = list(
                    values = parameter_grids$beta_m,
                    relative = beta_m_profile,
                    interval = beta_m_interval
                )
            ),
            reparameterised_profiles
        ),
        profile_summary = profile_summary,
        output_files = output_files
    )
}
